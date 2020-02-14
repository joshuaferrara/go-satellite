package satellite

import (
	"log"
	"math"
	"strconv"
	"strings"
)

const (
	twoPi   float64 = math.Pi * 2.0
	deg2rad float64 = math.Pi / 180.0
	xPDOTP  float64 = 1440.0 / (2.0 * math.Pi)
)

// LatLong holds latitude and Longitude in either degrees or radians
type LatLong struct {
	Latitude, Longitude float64
}

// Vector3 holds X, Y, Z position
type Vector3 struct {
	X, Y, Z float64
}

// LookAngles holds an azimuth, elevation and range
type LookAngles struct {
	Az, El, Rg float64
}

// ParseTLE parses a two line element dataset into a Satellite struct
func ParseTLE(line1, line2, gravconst string) (sat Satellite) {
	sat.Line1 = line1
	sat.Line2 = line2

	sat.Error = 0
	sat.whichconst = getGravConst(gravconst)

	// LINE 1 BEGIN
	sat.satnum = parseInt(strings.TrimSpace(line1[2:7]))
	sat.epochyr = parseInt(line1[18:20])
	sat.epochdays = parseFloat(line1[20:32])

	// These three can be negative / positive
	sat.ndot = parseFloat(strings.Replace(line1[33:43], " ", "", 2))
	sat.nddot = parseFloat(strings.Replace(line1[44:45]+"."+line1[45:50]+"e"+line1[50:52], " ", "", 2))
	sat.bstar = parseFloat(strings.Replace(line1[53:54]+"."+line1[54:59]+"e"+line1[59:61], " ", "", 2))
	// LINE 1 END

	// LINE 2 BEGIN
	sat.inclo = parseFloat(strings.Replace(line2[8:16], " ", "", 2))
	sat.nodeo = parseFloat(strings.Replace(line2[17:25], " ", "", 2))
	sat.ecco = parseFloat("." + line2[26:33])
	sat.argpo = parseFloat(strings.Replace(line2[34:42], " ", "", 2))
	sat.mo = parseFloat(strings.Replace(line2[43:51], " ", "", 2))
	sat.no = parseFloat(strings.Replace(line2[52:63], " ", "", 2))
	// LINE 2 END
	return
}

// TLEToSat converts a two line element data set into a Satellite struct and runs sgp4init
func TLEToSat(line1, line2 string, gravconst string) Satellite {
	//sat := Satellite{Line1: line1, Line2: line2}
	sat := ParseTLE(line1, line2, gravconst)

	opsmode := "i"

	sat.no = sat.no / xPDOTP
	sat.ndot = sat.ndot / (xPDOTP * 1440.0)
	sat.nddot = sat.nddot / (xPDOTP * 1440.0 * 1440)

	sat.inclo = sat.inclo * deg2rad
	sat.nodeo = sat.nodeo * deg2rad
	sat.argpo = sat.argpo * deg2rad
	sat.mo = sat.mo * deg2rad

	var year int64
	if sat.epochyr < 57 {
		year = sat.epochyr + 2000
	} else {
		year = sat.epochyr + 1900
	}

	mon, day, hr, min, sec := days2mdhms(year, sat.epochdays)

	sat.jdsatepoch = JDay(int(year), int(mon), int(day), int(hr), int(min), int(sec))

	sgp4init(&opsmode, sat.jdsatepoch-2433281.5, &sat)

	return sat
}

// Parses a string into a float64 value.
func parseFloat(strIn string) (ret float64) {
	ret, err := strconv.ParseFloat(strIn, 64)
	if err != nil {
		log.Fatal(err)
	}
	return ret
}

// Parses a string into a int64 value.
func parseInt(strIn string) (ret int64) {
	ret, err := strconv.ParseInt(strIn, 10, 0)
	if err != nil {
		log.Fatal(err)
	}
	return ret
}
