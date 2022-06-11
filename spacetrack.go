package satellite

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"net/http"
	"net/url"
	"strings"
	"time"

	"github.com/pkg/errors"
)

const authURL = "https://www.space-track.org/ajaxauth/login"
const baseurl = "https://www.space-track.org/basicspacedata/query/class"

type orderBy string

const asc orderBy = "asc"
const desc orderBy = "desc"

type operator string

const dq operator = ""
const lt operator = "<"

// Keeping these around in case the API interface expands to use them
// const gt operator = ">"
// const ne operator = "<>"
// const or operator = ","
// const inclusiveRange operator = "--"
// const null operator = "null-val"
// const like operator = "~~"
// const wildcard operator = "^"
// const now operator = "now"

type field string

const epoch field = "EPOCH"
const tle field = "TLE"
const noradCatID field = "NORAD_CAT_ID"

var ErrInvalidResponseCode = errors.New("Invalid response from spacetrack")
var ErrNotSingleSat = errors.New("not a single satellite returned")

// Spacetrack contains an initialised API interface to space-track.org
type Spacetrack struct {
	username string
	password string
}

// NewSpacetrack creates an initialised API interface to space-track.org
// https://space-track.org allows you to create a free account, however be
// aware there *are* API throttles (see https://www.space-track.org/documentation#/api )
func NewSpacetrack(username, password string) *Spacetrack {
	return &Spacetrack{
		username: username,
		password: password,
	}
}

// GetTLE generates a Satellite from the most recent TLE from space-track.org before the given time
func (s *Spacetrack) GetTLE(catid uint64, ts time.Time, gravConst Gravity) (Satellite, error) {
	zero := Satellite{}
	args := spacetrackArgs{
		base:         baseurl,
		class:        tle,
		orderByField: epoch,
		orderByDir:   desc,
		limit:        1,
		filters: filters{{
			filterType: noradCatID,
			operator:   dq,
			value:      fmt.Sprint(catid),
		}, {
			filterType: epoch,
			operator:   lt,
			value:      ts.Format(time.RFC3339),
		}},
	}
	q := buildQuery(args)

	vals := url.Values{}
	vals.Add("identity", s.username)
	vals.Add("password", s.password)
	vals.Add("query", q)

	client := &http.Client{}

	resp, err := client.PostForm(authURL, vals)
	if err != nil {
		return zero, err
	}
	defer resp.Body.Close()
	if resp.StatusCode < 200 || resp.StatusCode >= 300 {
		return zero, errors.Wrap(ErrInvalidResponseCode, resp.Status)
	}
	respData, err := ioutil.ReadAll(resp.Body)
	if err != nil {
		return zero, err
	}

	var sats []Satellite
	if err := json.Unmarshal(respData, &sats); err != nil {
		return zero, err
	}
	if len(sats) != 1 {
		return zero, errors.Wrap(ErrNotSingleSat, fmt.Sprint(sats))

	}
	sat := TLEToSat(sats[0].Line1, sats[0].Line2, gravConst)
	return sat, nil
}

type spacetrackArgs struct {
	base         string
	class        field
	orderByField field
	orderByDir   orderBy
	limit        uint64
	filters      filters
}

func (s spacetrackArgs) orderBy() string {
	if s.orderByField == "" {
		return ""
	}
	if s.orderByDir == "" {
		s.orderByDir = asc
	}
	return fmt.Sprint("orderby/", s.orderByField, " ", s.orderByDir)
}

func (s spacetrackArgs) limitString() string {
	if s.limit == 0 {
		return ""
	}
	return fmt.Sprint("limit/", s.limit)
}

type filters []filter
type filter struct {
	filterType field
	operator   operator
	value      string
}

func (f filters) render() string {
	results := []string{}
	for _, f := range f {
		results = append(results, f.render())
	}
	return strings.Join(results, "/")
}

func (f filter) render() string {
	return fmt.Sprint(f.filterType, "/", f.operator, f.value)
}

func buildQuery(a spacetrackArgs) string {
	fields := []string{
		a.base,
		string(a.class),
	}
	optionalFields := []string{
		a.filters.render(),
		a.orderBy(),
		a.limitString(),
		"emptyresult/show",
	}

	for _, f := range optionalFields {
		if f == "" {
			continue
		}
		fields = append(fields, f)
	}

	return strings.Join(fields, "/")

}
