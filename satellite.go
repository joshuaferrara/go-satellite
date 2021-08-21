package satellite

// Struct for holding satellite information during and before propagation
type Satellite struct {
	Line1 string `json:"TLE_LINE1"`
	Line2 string `json:"TLE_LINE2"`

	satnum int64

	Error      int64
	ErrorStr   string
	whichconst GravConst

	epochyr    int64
	epochdays  float64
	jdsatepoch float64

	ndot  float64
	nddot float64
	bstar float64
	inclo float64
	nodeo float64
	ecco  float64
	argpo float64
	mo    float64
	no    float64
	alta  float64
	altp  float64

	method        string
	operationmode string
	init          string

	gsto    float64
	isimp   float64
	con41   float64
	cc5     float64
	d4      float64
	argpdot float64
	t       float64
	t4cof   float64
	x7thm1  float64
	xlcof   float64
	cc1     float64
	d2      float64
	delmo   float64
	omgcof  float64
	t2cof   float64
	t5cof   float64
	mdot    float64
	xmcof   float64
	aycof   float64
	cc4     float64
	d3      float64
	eta     float64
	sinmao  float64
	t3cof   float64
	x1mth2  float64
	nodedot float64
	nodecf  float64

	irez  float64
	d3210 float64
	d4422 float64
	d5421 float64
	del1  float64
	didt  float64
	domdt float64
	peo   float64
	pinco float64
	se3   float64
	sgh4  float64
	si2   float64
	sl3   float64
	xfact float64
	xgh4  float64
	xi2   float64
	xl3   float64
	zmol  float64
	xli   float64
	d2201 float64
	d3222 float64
	d5220 float64
	d5433 float64
	del2  float64
	dmdt  float64
	e3    float64
	pgho  float64
	plo   float64
	sgh2  float64
	sh2   float64
	si3   float64
	sl4   float64
	xgh2  float64
	xh2   float64
	xi3   float64
	xl4   float64
	zmos  float64
	xni   float64
	d2211 float64
	d4410 float64
	d5232 float64
	dedt  float64
	del3  float64
	dnodt float64
	ee2   float64
	pho   float64
	se2   float64
	sgh3  float64
	sh3   float64
	sl2   float64
	xgh3  float64
	xh3   float64
	xl2   float64
	xlamo float64
	atime float64
}
