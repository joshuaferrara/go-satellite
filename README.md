# satellite
    import "github.com/joshuaferrara/go-satellite"

## Intro

[![Go](https://github.com/joshuaferrara/go-satellite/actions/workflows/go.yml/badge.svg?branch=master)](https://github.com/joshuaferrara/go-satellite/actions/workflows/go.yml) [![GoDoc](https://godoc.org/github.com/joshuaferrara/go-satellite?status.svg)](https://godoc.org/github.com/joshuaferrara/go-satellite)

I decided to port the SGP4 library to GoLang as one of my first projects with the language. I've included a test suite to ensure accuracy.

## Usage

#### Constants

```go
const DEG2RAD float64 = math.Pi / 180.0
```

```go
const RAD2DEG float64 = 180.0 / math.Pi
```

```go
const TWOPI float64 = math.Pi * 2.0
```

```go
const XPDOTP float64 = 1440.0 / (2.0 * math.Pi)
```

#### func  ECIToLLA

```go
func ECIToLLA(eciCoords Vector3, gmst float64) (altitude, velocity float64, ret LatLong)
```
Convert Earth Centered Inertial coordinated into equivalent latitude, longitude,
altitude and velocity. Reference: http://celestrak.com/columns/v02n03/

#### func  GSTimeFromDate

```go
func GSTimeFromDate(year, mon, day, hr, min, sec int) float64
```
Calc GST given year, month, day, hour, minute and second

#### func  JDay

```go
func JDay(year, mon, day, hr, min, sec int) float64
```
Calc julian date given year, month, day, hour, minute and second the julian date
is defined by each elapsed day since noon, jan 1, 4713 bc.

#### func  Propagate

```go
func Propagate(sat Satellite, year int, month int, day, hours, minutes, seconds int) (position, velocity Vector3)
```
Calculates position and velocity vectors for given time

#### func  ThetaG_JD

```go
func ThetaG_JD(jday float64) (ret float64)
```
Calculate GMST from Julian date. Reference: The 1992 Astronomical Almanac, page
B6.


#### type LatLong

```go
type LatLong struct {
	Latitude, Longitude float64
}
```

Holds latitude and Longitude in either degrees or radians

#### func  LatLongDeg

```go
func LatLongDeg(rad LatLong) (deg LatLong)
```
Convert LatLong in radians to LatLong in degrees

#### type LookAngles

```go
type LookAngles struct {
	Az, El, Rg float64
}
```

Holds an azimuth, elevation and range

#### func  ECIToLookAngles

```go
func ECIToLookAngles(eciSat Vector3, obsCoords LatLong, obsAlt, jday float64) (lookAngles LookAngles)
```
Calculate look angles for given satellite position and observer position obsAlt
in km Reference: http://celestrak.com/columns/v02n02/

#### type Satellite

```go
type Satellite struct {
	Line1 string
	Line2 string
}
```

Struct for holding satellite information during and before propagation

#### func  ParseTLE

```go
func ParseTLE(line1, line2 string, gravConst Gravity) (sat Satellite)
```
Parses a two line element dataset into a Satellite struct

#### func  TLEToSat

```go
func TLEToSat(line1, line2 string, gravConst Gravity) Satellite
```
Converts a two line element data set into a Satellite struct and runs sgp4init

#### type Vector3

```go
type Vector3 struct {
	X, Y, Z float64
}
```

Holds X, Y, Z position

#### func  ECIToECEF

```go
func ECIToECEF(eciCoords Vector3, gmst float64) (ecfCoords Vector3)
```
Convert Earth Centered Intertial coordinates into Earth Cenetered Earth Final
coordinates Reference: http://ccar.colorado.edu/ASEN5070/handouts/coordsys.doc

#### func  LLAToECI

```go
func LLAToECI(obsCoords LatLong, alt, jday float64) (eciObs Vector3)
```
Convert latitude, longitude and altitude into equivalent Earth Centered
Intertial coordinates Reference: The 1992 Astronomical Almanac, page K11.

#### func  NewSpacetrack
```go
func NewSpacetrack(username, password string) *Spacetrack
```
Initialise a spacetrack API for fetching TLEs

#### func  Spacetrack.GetTLE()
```go
func (s *Spacetrack) GetTLE(catid uint64, ts time.Time, gravConst Gravity) (Satellite, error)
```
Get an initialized Satellite based on the latest TLE before the given time.
