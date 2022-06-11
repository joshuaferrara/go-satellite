package satellite

import (
	"log"
	"math"
)

// Holds variables that are dependent upon selected gravity model
type GravConst struct {
	mu, radiusearthkm, xke, tumin, j2, j3, j4, j3oj2 float64
}

type Gravity string

const (
	GravityWGS72Old Gravity = "wgs72old"
	GravityWGS72    Gravity = "wgs72"
	GravityWGS84    Gravity = "wgs84"
)

// Returns a GravConst with correct information on requested model provided through the name parameter
func getGravConst(name Gravity) (grav GravConst) {
	switch name {
	case GravityWGS72Old:
		grav.mu = 398600.79964
		grav.radiusearthkm = 6378.135
		grav.xke = 0.0743669161
		grav.tumin = 1.0 / grav.xke
		grav.j2 = 0.001082616
		grav.j3 = -0.00000253881
		grav.j4 = -0.00000165597
		grav.j3oj2 = grav.j3 / grav.j2
	case GravityWGS72:
		grav.mu = 398600.8
		grav.radiusearthkm = 6378.135
		grav.xke = 60.0 / math.Sqrt(grav.radiusearthkm*grav.radiusearthkm*grav.radiusearthkm/grav.mu)
		grav.tumin = 1.0 / grav.xke
		grav.j2 = 0.001082616
		grav.j3 = -0.00000253881
		grav.j4 = -0.00000165597
		grav.j3oj2 = grav.j3 / grav.j2
	case GravityWGS84:
		grav.mu = 398600.5
		grav.radiusearthkm = 6378.137
		grav.xke = 60.0 / math.Sqrt(grav.radiusearthkm*grav.radiusearthkm*grav.radiusearthkm/grav.mu)
		grav.tumin = 1.0 / grav.xke
		grav.j2 = 0.00108262998905
		grav.j3 = -0.00000253215306
		grav.j4 = -0.00000161098761
		grav.j3oj2 = grav.j3 / grav.j2
	default:
		log.Fatal(name, "is not a valid gravity model")
	}

	return
}

// Not the movie
