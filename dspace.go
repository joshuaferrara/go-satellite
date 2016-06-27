package satellite

import (
	"math"
)

// A struct returned from the dsinit function
type DeepSpaceInitResult struct {
	em, argpm, inclm, mm, nm, nodem, irez, atime, d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433, dedt, didt, dmdt, dndt, dnodt, domdt, del1, del2, del3, xfact, xlamo, xli, xni float64
}

// A struct returned from the dspace function
type DeepSpaceResult struct {
	atime, em, argpm, inclm, xli, mm, xni, nodem, dndt, nm float64
}

// this procedure provides deep space contributions to mean motion dot due to geopotential resonance with half day and one day orbits.
func dsinit(whichconst GravConst, cosim, emsq, argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4, ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, t, tc, gsto, mo, mdot, no, nodeo, nodedot, xpidot, z1, z3, z11, z13, z21, z23, z31, z33, ecco, eccsq, em, argpm, inclm, mm, nm, nodem, irez, atime, d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433, dedt, didt, dmdt, dnodt, domdt, del1, del2, del3, xfact, xlamo, xli, xni float64) (res DeepSpaceInitResult) {

	var f220, f221, f311, f321, f322, f330, f441, f442, f522, f523, f542, f543, g200, g201, g211, g300, g310, g322, g410, g422, g520, g521, g532, g533, sini2, temp, temp1, theta, xno2, ainv2, aonv, cosisq, eoc float64

	q22 := 1.7891679e-6
	q31 := 2.1460748e-6
	q33 := 2.2123015e-7
	root22 := 1.7891679e-6
	root44 := 7.3636953e-9
	root54 := 2.1765803e-9
	rptim := 4.37526908801129966e-3
	root32 := 3.7393792e-7
	root52 := 1.1428639e-7
	x2o3 := 2.0 / 3.0
	znl := 1.5835218e-4
	zns := 1.19459e-5

	xke := whichconst.xke

	irez = 0
	if 0.0034906585 < nm && nm < 0.0052359877 {
		irez = 1
	}
	if 8.26e-3 <= nm && nm <= 9.24e-3 && em >= 0.5 {
		irez = 2
	}

	ses := ss1 * zns * ss5
	sis := ss2 * zns * (sz11 + sz13)
	sls := -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0*emsq)
	sghs := ss4 * zns * (sz31 + sz33 - 6.0)
	shs := -zns * ss2 * (sz21 + sz23)

	if inclm < 5.2359877e-2 || inclm > math.Pi-5.2359877e-2 {
		shs = 0.0
	}
	if sinim != 0.0 {
		shs = shs / sinim
	}
	sgs := sghs - cosim*shs

	dedt = ses + s1*znl*s5
	didt = sis + s2*znl*(z11+z13)
	dmdt = sls - znl*s3*(z1+z3-14.0-6.0*emsq)
	sghl := s4 * znl * (z31 + z33 - 6.0)
	shll := -znl * s2 * (z21 + z23)

	if inclm < 5.2359877e-2 || inclm > math.Pi-5.2359877e-2 {
		shll = 0.0
	}
	domdt = sgs + sghl
	dnodt = shs

	if sinim != 0.0 {
		domdt = domdt - cosim/sinim*shll
		dnodt = dnodt + shll/sinim
	}

	dndt := 0.0
	theta = math.Mod(gsto+tc*rptim, TWOPI)
	em = em + dedt*t
	inclm = inclm + didt*t
	argpm = argpm + domdt*t
	nodem = nodem + dnodt*t
	mm = mm + dmdt*t

	if irez != 0 {
		aonv = math.Pow(nm/xke, x2o3)
		if irez == 2 {
			cosisq = cosim * cosim
			emo := em
			em = ecco
			emsqo := emsq
			emsq = eccsq
			eoc = em * emsq
			g201 = -0.306 - (em-0.64)*0.440

			if em <= 0.65 {
				g211 = 3.616 - 13.2470*em + 16.2900*emsq
				g310 = -19.302 + 117.3900*em - 228.4190*emsq + 156.5910*eoc
				g322 = -18.9068 + 109.7927*em - 214.6334*emsq + 146.5816*eoc
				g410 = -41.122 + 242.6940*em - 471.0940*emsq + 313.9530*eoc
				g422 = -146.407 + 841.8800*em - 1629.014*emsq + 1083.4350*eoc
				g520 = -532.114 + 3017.977*em - 5740.032*emsq + 3708.2760*eoc
			} else {
				g211 = -72.099 + 331.819*em - 508.738*emsq + 266.724*eoc
				g310 = -346.844 + 1582.851*em - 2415.925*emsq + 1246.113*eoc
				g322 = -342.585 + 1554.908*em - 2366.899*emsq + 1215.972*eoc
				g410 = -1052.797 + 4758.686*em - 7193.992*emsq + 3651.957*eoc
				g422 = -3581.690 + 16178.110*em - 24462.770*emsq + 12422.520*eoc
				if em > 0.715 {
					g520 = -5149.66 + 29936.92*em - 54087.36*emsq + 31324.56*eoc
				} else {
					g520 = 1464.74 - 4664.75*em + 3763.64*emsq
				}
			}

			if em < 0.7 {
				g533 = -919.22770 + 4988.6100*em - 9064.7700*emsq + 5542.21*eoc
				g521 = -822.71072 + 4568.6173*em - 8491.4146*emsq + 5337.524*eoc
				g532 = -853.66600 + 4690.2500*em - 8624.7700*emsq + 5341.4*eoc
			} else {
				g533 = -37995.780 + 161616.52*em - 229838.20*emsq + 109377.94*eoc
				g521 = -51752.104 + 218913.95*em - 309468.16*emsq + 146349.42*eoc
				g532 = -40023.880 + 170470.89*em - 242699.48*emsq + 115605.82*eoc
			}

			sini2 = sinim * sinim
			f220 = 0.75 * (1.0 + 2.0*cosim + cosisq)
			f221 = 1.5 * sini2
			f321 = 1.875 * sinim * (1.0 - 2.0*cosim - 3.0*cosisq)
			f322 = -1.875 * sinim * (1.0 + 2.0*cosim - 3.0*cosisq)
			f441 = 35.0 * sini2 * f220
			f442 = 39.3750 * sini2 * sini2
			f522 = 9.84375 * sinim * (sini2*(1.0-2.0*cosim-5.0*cosisq) + 0.33333333*(-2.0+4.0*cosim+6.0*cosisq))
			f523 = sinim * (4.92187512*sini2*(-2.0-4.0*cosim+10.0*cosisq) + 6.56250012*(1.0+2.0*cosim-3.0*cosisq))
			f542 = 29.53125 * sinim * (2.0 - 8.0*cosim + cosisq*(-12.0+8.0*cosim+10.0*cosisq))
			f543 = 29.53125 * sinim * (-2.0 - 8.0*cosim + cosisq*(12.0+8.0*cosim-10.0*cosisq))
			xno2 = nm * nm
			ainv2 = aonv * aonv
			temp1 = 3.0 * xno2 * ainv2
			temp = temp1 * root22
			d2201 = temp * f220 * g201
			d2211 = temp * f221 * g211
			temp1 = temp1 * aonv
			temp = temp1 * root32
			d3210 = temp * f321 * g310
			d3222 = temp * f322 * g322
			temp1 = temp1 * aonv
			temp = 2.0 * temp1 * root44
			d4410 = temp * f441 * g410
			d4422 = temp * f442 * g422
			temp1 = temp1 * aonv
			temp = temp1 * root52
			d5220 = temp * f522 * g520
			d5232 = temp * f523 * g532
			temp = 2.0 * temp1 * root54
			d5421 = temp * f542 * g521
			d5433 = temp * f543 * g533
			xlamo = math.Mod(mo+nodeo+nodeo-theta-theta, TWOPI)
			xfact = mdot + dmdt + 2.0*(nodedot+dnodt-rptim) - no
			em = emo
			emsq = emsqo
		}
		if irez == 1 {
			g200 = 1.0 + emsq*(-2.5+0.8125*emsq)
			g310 = 1.0 + 2.0*emsq
			g300 = 1.0 + emsq*(-6.0+6.60937*emsq)
			f220 = 0.75 * (1.0 + cosim) * (1.0 + cosim)
			f311 = 0.9375*sinim*sinim*(1.0+3.0*cosim) - 0.75*(1.0+cosim)
			f330 = 1.0 + cosim
			f330 = 1.875 * f330 * f330 * f330
			del1 = 3.0 * nm * nm * aonv * aonv
			del2 = 2.0 * del1 * f220 * g200 * q22
			del3 = 3.0 * del1 * f330 * g300 * q33 * aonv
			del1 = del1 * f311 * g310 * q31 * aonv
			xlamo = math.Mod(mo+nodeo+argpo-theta, TWOPI)
			xfact = mdot + xpidot - rptim + dmdt + domdt + dnodt - no
		}
		xli = xlamo
		xni = no
		atime = 0.0
		nm = no + dndt
	}

	res.em = em
	res.argpm = argpm
	res.inclm = inclm
	res.mm = mm
	res.nm = nm
	res.nodem = nodem
	res.irez = irez
	res.atime = atime
	res.d2201 = d2201
	res.d2211 = d2211
	res.d3210 = d3210
	res.d3222 = d3222
	res.d4410 = d4410
	res.d4422 = d4422
	res.d5220 = d5220
	res.d5232 = d5232
	res.d5421 = d5421
	res.d5433 = d5433
	res.dedt = dedt
	res.didt = didt
	res.dmdt = dmdt
	res.dndt = dndt
	res.dnodt = dnodt
	res.domdt = domdt
	res.del1 = del1
	res.del2 = del2
	res.del3 = del3
	res.xfact = xfact
	res.xlamo = xlamo
	res.xli = xli
	res.xni = xni

	return
}

// this procedure provides deep space contributions to mean elements for perturbing third body. these effects have been averaged over one revolution of the sun and moon. for earth resonance effects, the effects have been averaged over no revolutions of the satellite. (mean motion)
func dspace(irez, d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433, dedt, del1, del2, del3, didt, dmdt, dnodt, domdt, argpo, argpdot, t, tc, gsto, xfact, xlamo, no, atime, em, argpm, inclm, xli, mm, xni, nodem, nm float64) (result DeepSpaceResult) {
	var delt, ft, theta, x2li, x2omi, xl, xldot, xnddt, xndt, xomi float64

	fasx2 := 0.13130908
	fasx4 := 2.8843198
	fasx6 := 0.37448087
	g22 := 5.7686396
	g32 := 0.95240898
	g44 := 1.8014998
	g52 := 1.0508330
	g54 := 4.4108898
	rptim := 4.37526908801129966e-3
	stepp := 720.0
	stepn := -720.0
	step2 := 259200.0

	dndt := 0.0
	theta = math.Mod((gsto + tc*rptim), TWOPI)
	em = em + dedt*t

	inclm = inclm + didt*t
	argpm = argpm + domdt*t
	nodem = nodem + dnodt*t
	mm = mm + dmdt*t

	ft = 0.0

	if irez != 0.0 {
		if atime == 0.0 || t*atime <= 0.0 || math.Abs(t) < math.Abs(atime) {
			atime = 0.0
			xni = no
			xli = xlamo
		}

		if t > 0.0 {
			delt = stepp
		} else {
			delt = stepn
		}

		iretn := 381
		for iretn == 381 {
			if irez != 2 {
				xndt = del1*math.Sin(xli-fasx2) + del2*math.Sin(2.0*(xli-fasx4)) + del3*math.Sin(3.0*(xli-fasx6))
				xldot = xni + xfact
				xnddt = del1*math.Cos(xli-fasx2) + 2.0*del2*math.Cos(2.0*(xli-fasx4)) + 3.0*del3*math.Cos(3.0*(xli-fasx6))
				xnddt = xnddt * xldot
			} else {
				xomi = argpo + argpdot*atime
				x2omi = xomi + xomi
				x2li = xli + xli
				xndt = (d2201*math.Sin(x2omi+xli-g22) + d2211*math.Sin(xli-g22) + d3210*math.Sin(xomi+xli-g32) + d3222*math.Sin(-xomi+xli-g32) + d4410*math.Sin(x2omi+x2li-g44) + d4422*math.Sin(x2li-g44) + d5220*math.Sin(xomi+xli-g52) + d5232*math.Sin(-xomi+xli-g52) + d5421*math.Sin(xomi+x2li-g54) + d5433*math.Sin(-xomi+x2li-g54))
				xldot = xni + xfact
				xnddt = (d2201*math.Cos(x2omi+xli-g22) + d2211*math.Cos(xli-g22) + d3210*math.Cos(xomi+xli-g32) + d3222*math.Cos(-xomi+xli-g32) + d5220*math.Cos(xomi+xli-g52) + d5232*math.Cos(-xomi+xli-g52) + 2.0*(d4410*math.Cos(x2omi+x2li-g44)+d4422*math.Cos(x2li-g44)+d5421*math.Cos(xomi+x2li-g54)+d5433*math.Cos(-xomi+x2li-g54)))
				xnddt = xnddt * xldot
			}

			if math.Abs(t-atime) >= stepp {
				iretn = 381
			} else {
				ft = t - atime
				iretn = 0
			}

			if iretn == 381 {
				xli = xli + xldot*delt + xndt*step2
				xni = xni + xndt*delt + xnddt*step2
				atime = atime + delt
			}
		}

		nm = xni + xndt*ft + xnddt*ft*ft*0.5
		xl = xli + xldot*ft + xndt*ft*ft*0.5
		if irez != 1 {
			mm = xl - 2.0*nodem + 2.0*theta
			dndt = nm - no
		} else {
			mm = xl - nodem - argpm + theta
			dndt = nm - no
		}

		nm = no + dndt
	}

	result.atime = atime
	result.em = em
	result.argpm = argpm
	result.inclm = inclm
	result.xli = xli
	result.mm = mm
	result.xni = xni
	result.nodem = nodem
	result.dndt = dndt
	result.nm = nm

	return
}

// A struct returned from the dpper function
type DpperResult struct {
	ep, inclp, nodep, argpp, mp float64
}

// this procedure provides deep space long period periodic contributions to the mean elements. by design, these periodics are zero at epoch. this used to be dscom which included initialization, but it's really a recurring function.
func dpper(satrec *Satellite, inclo float64, init string, ep, inclp, nodep, argpp, mp float64, opsmode string) (result DpperResult) {
	e3 := satrec.e3
	ee2 := satrec.ee2
	peo := satrec.peo
	pgho := satrec.pgho
	pho := satrec.pho
	pinco := satrec.pinco
	plo := satrec.plo
	se2 := satrec.se2
	se3 := satrec.se3
	sgh2 := satrec.sgh2
	sgh3 := satrec.sgh3
	sgh4 := satrec.sgh4
	sh2 := satrec.sh2
	sh3 := satrec.sh3
	si2 := satrec.si2
	si3 := satrec.si3
	sl2 := satrec.sl2
	sl3 := satrec.sl3
	sl4 := satrec.sl4
	t := satrec.t
	xgh2 := satrec.xgh2
	xgh3 := satrec.xgh3
	xgh4 := satrec.xgh4
	xh2 := satrec.xh2
	xh3 := satrec.xh3
	xi2 := satrec.xi2
	xi3 := satrec.xi3
	xl2 := satrec.xl2
	xl3 := satrec.xl3
	xl4 := satrec.xl4
	zmol := satrec.zmol
	zmos := satrec.zmos

	zns := 1.19459e-5
	zes := 0.01675
	znl := 1.5835218e-4
	zel := 0.05490

	zm := zmos + zns*t

	if init == "y" {
		zm = zmos
	}

	zf := zm + 2.0*zes*math.Sin(zm)
	sinzf := math.Sin(zf)
	f2 := 0.5*sinzf*sinzf - 0.25
	f3 := -0.5 * sinzf * math.Cos(zf)
	ses := se2*f2 + se3*f3
	sis := si2*f2 + si3*f3
	sls := sl2*f2 + sl3*f3 + sl4*sinzf
	sghs := sgh2*f2 + sgh3*f3 + sgh4*sinzf
	shs := sh2*f2 + sh3*f3

	zm = zmol + znl*t
	if init == "y" {
		zm = zmol
	}

	zf = zm + 2.0*zel*math.Sin(zm)
	sinzf = math.Sin(zf)
	f2 = 0.5*sinzf*sinzf - 0.25
	f3 = -0.5 * sinzf * math.Cos(zf)

	sel := ee2*f2 + e3*f3
	sil := xi2*f2 + xi3*f3
	sll := xl2*f2 + xl3*f3 + xl4*sinzf
	sghl := xgh2*f2 + xgh3*f3 + xgh4*sinzf
	shll := xh2*f2 + xh3*f3
	pe := ses + sel
	pinc := sis + sil
	pl := sls + sll
	pgh := sghs + sghl
	ph := shs + shll

	if init == "n" {
		pe = pe - peo
		pinc = pinc - pinco
		pl = pl - plo
		pgh = pgh - pgho
		ph = ph - pho
		inclp = inclp + pinc
		ep = ep + pe
		sinip := math.Sin(inclp)
		cosip := math.Cos(inclp)

		if inclp >= 0.2 {
			ph /= sinip
			pgh -= cosip * ph
			argpp += pgh
			nodep += ph
			mp += pl
		} else {
			sinop := math.Sin(nodep)
			cosop := math.Cos(nodep)
			alfdp := sinip * sinop
			betdp := sinip * cosop
			dalf := ph*cosop + pinc*cosip*sinop
			dbet := -ph*sinop + pinc*cosip*cosop

			alfdp = alfdp + dalf
			betdp = betdp + dbet
			nodep = math.Mod(nodep, TWOPI)
			if nodep < 0.0 && opsmode == "a" {
				nodep = nodep + TWOPI
			}
			xls := mp + argpp + pl + pgh + (cosip-pinc*sinip)*nodep
			xnoh := nodep
			nodep = math.Atan2(alfdp, betdp)
			if nodep < 0.0 && opsmode == "a" {
				nodep = nodep + TWOPI
			}
			if math.Abs(xnoh-nodep) > math.Pi {
				if nodep < xnoh {
					nodep = nodep + TWOPI
				} else {
					nodep = nodep - TWOPI
				}
			}
			mp += pl
			argpp = xls - mp - cosip*nodep
		}
	}

	result.ep = ep
	result.inclp = inclp
	result.nodep = nodep
	result.argpp = argpp
	result.mp = mp

	return
}

// A struct returned from the dscom function
type DSComResults struct {
	snodm, cnodm, sinim, cosim, sinomm, cosomm, day, e3, ee2, em, emsq, gam, peo, pgho, pho, pinco, plo, rtemsq, se2, se3, sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4, s1, s2, s3, s4, s5, s6, s7, ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3, sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33, xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4, nm, z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33, zmol, zmos float64
}

// this procedure provides deep space common items used by both the secular and periodics subroutines. input is provided as shown. this routine used to be called dpper, but the functions inside weren't well organized.
func dscom(epoch, ep, argpp, tc, inclp, nodep, np, e3, ee2, peo, pgho, pho, pinco, plo, se2, se3, sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4, xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4, zmol, zmos float64) (res DSComResults) {
	var a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, betasq, cc, ctem, stem, x1, x2, x3, x4, x5, x6, x7, x8, xnodce, xnoi, zcosg, zsing, zcosgl, zsingl, zcosh, zsinh, zcoshl, zsinhl, zcosi, zsini, zcosil, zsinil, zx, zy, ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3, sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33, s1, s2, s3, s4, s5, s6, s7, z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33 float64

	zes := 0.01675
	zel := 0.05490
	c1ss := 2.9864797e-6
	c1l := 4.7968065e-7
	zsinis := 0.39785416
	zcosis := 0.91744867
	zcosgs := 0.1945905
	zsings := -0.98088458

	nm := np
	em := ep
	snodm := math.Sin(nodep)
	cnodm := math.Cos(nodep)
	sinomm := math.Sin(argpp)
	cosomm := math.Cos(argpp)
	sinim := math.Sin(inclp)
	cosim := math.Cos(inclp)
	emsq := em * em
	betasq = 1.0 - emsq
	rtemsq := math.Sqrt(betasq)

	peo = 0.0
	pinco = 0.0
	plo = 0.0
	pgho = 0.0
	pho = 0.0
	day := epoch + 18261.5 + tc/1440.0
	xnodce = math.Mod(4.5236020-9.2422029e-4*day, TWOPI)
	stem = math.Sin(xnodce)
	ctem = math.Cos(xnodce)
	zcosil = 0.91375164 - 0.03568096*ctem
	zsinil = math.Sqrt(1.0 - zcosil*zcosil)
	zsinhl = 0.089683511 * stem / zsinil
	zcoshl = math.Sqrt(1.0 - zsinhl*zsinhl)
	gam := 5.8351514 + 0.0019443680*day
	zx = 0.39785416 * stem / zsinil
	zy = zcoshl*ctem + 0.91744867*zsinhl*stem
	zx = math.Atan2(zx, zy)
	zx = gam + zx - xnodce
	zcosgl = math.Cos(zx)
	zsingl = math.Sin(zx)

	zcosg = zcosgs
	zsing = zsings
	zcosi = zcosis
	zsini = zsinis
	zcosh = cnodm
	zsinh = snodm
	cc = c1ss
	xnoi = 1.0 / nm

	for lsflg := 1; lsflg <= 2; lsflg++ {
		a1 = zcosg*zcosh + zsing*zcosi*zsinh
		a3 = -zsing*zcosh + zcosg*zcosi*zsinh
		a7 = -zcosg*zsinh + zsing*zcosi*zcosh
		a8 = zsing * zsini
		a9 = zsing*zsinh + zcosg*zcosi*zcosh
		a10 = zcosg * zsini
		a2 = cosim*a7 + sinim*a8
		a4 = cosim*a9 + sinim*a10
		a5 = -sinim*a7 + cosim*a8
		a6 = -sinim*a9 + cosim*a10

		x1 = a1*cosomm + a2*sinomm
		x2 = a3*cosomm + a4*sinomm
		x3 = -a1*sinomm + a2*cosomm
		x4 = -a3*sinomm + a4*cosomm
		x5 = a5 * sinomm
		x6 = a6 * sinomm
		x7 = a5 * cosomm
		x8 = a6 * cosomm

		z31 = 12.0*x1*x1 - 3.0*x3*x3
		z32 = 24.0*x1*x2 - 6.0*x3*x4
		z33 = 12.0*x2*x2 - 3.0*x4*x4
		z1 = 3.0*(a1*a1+a2*a2) + z31*emsq
		z2 = 6.0*(a1*a3+a2*a4) + z32*emsq
		z3 = 3.0*(a3*a3+a4*a4) + z33*emsq
		z11 = -6.0*a1*a5 + emsq*(-24.0*x1*x7-6.0*x3*x5)
		z12 = -6.0*(a1*a6+a3*a5) + emsq*(-24.0*(x2*x7+x1*x8)-6.0*(x3*x6+x4*x5))
		z13 = -6.0*a3*a6 + emsq*(-24.0*x2*x8-6.0*x4*x6)
		z21 = 6.0*a2*a5 + emsq*(24.0*x1*x5-6.0*x3*x7)
		z22 = 6.0*(a4*a5+a2*a6) + emsq*(24.0*(x2*x5+x1*x6)-6.0*(x4*x7+x3*x8))
		z23 = 6.0*a4*a6 + emsq*(24.0*x2*x6-6.0*x4*x8)
		z1 = z1 + z1 + betasq*z31
		z2 = z2 + z2 + betasq*z32
		z3 = z3 + z3 + betasq*z33
		s3 = cc * xnoi
		s2 = -0.5 * s3 / rtemsq
		s4 = s3 * rtemsq
		s1 = -15.0 * em * s4
		s5 = x1*x3 + x2*x4
		s6 = x2*x3 + x1*x4
		s7 = x2*x4 - x1*x3

		if lsflg == 1 {
			ss1 = s1
			ss2 = s2
			ss3 = s3
			ss4 = s4
			ss5 = s5
			ss6 = s6
			ss7 = s7
			sz1 = z1
			sz2 = z2
			sz3 = z3
			sz11 = z11
			sz12 = z12
			sz13 = z13
			sz21 = z21
			sz22 = z22
			sz23 = z23
			sz31 = z31
			sz32 = z32
			sz33 = z33
			zcosg = zcosgl
			zsing = zsingl
			zcosi = zcosil
			zsini = zsinil
			zcosh = zcoshl*cnodm + zsinhl*snodm
			zsinh = snodm*zcoshl - cnodm*zsinhl
			cc = c1l
		}
	}

	zmol = math.Mod(4.7199672+0.22997150*day-gam, TWOPI)
	zmos = math.Mod(6.2565837+0.017201977*day, TWOPI)

	se2 = 2.0 * ss1 * ss6
	se3 = 2.0 * ss1 * ss7
	si2 = 2.0 * ss2 * sz12
	si3 = 2.0 * ss2 * (sz13 - sz11)
	sl2 = -2.0 * ss3 * sz2
	sl3 = -2.0 * ss3 * (sz3 - sz1)
	sl4 = -2.0 * ss3 * (-21.0 - 9.0*emsq) * zes
	sgh2 = 2.0 * ss4 * sz32
	sgh3 = 2.0 * ss4 * (sz33 - sz31)
	sgh4 = -18.0 * ss4 * zes
	sh2 = -2.0 * ss2 * sz22
	sh3 = -2.0 * ss2 * (sz23 - sz21)

	ee2 = 2.0 * s1 * s6
	e3 = 2.0 * s1 * s7
	xi2 = 2.0 * s2 * z12
	xi3 = 2.0 * s2 * (z13 - z11)
	xl2 = -2.0 * s3 * z2
	xl3 = -2.0 * s3 * (z3 - z1)
	xl4 = -2.0 * s3 * (-21.0 - 9.0*emsq) * zel
	xgh2 = 2.0 * s4 * z32
	xgh3 = 2.0 * s4 * (z33 - z31)
	xgh4 = -18.0 * s4 * zel
	xh2 = -2.0 * s2 * z22
	xh3 = -2.0 * s2 * (z23 - z21)

	res.snodm = snodm
	res.cnodm = cnodm
	res.sinim = sinim
	res.cosim = cosim
	res.sinomm = sinomm
	res.cosomm = cosomm
	res.day = day
	res.e3 = e3
	res.ee2 = ee2
	res.em = em
	res.emsq = emsq
	res.gam = gam
	res.peo = peo
	res.pgho = pgho
	res.pho = pho
	res.pinco = pinco
	res.plo = plo
	res.rtemsq = rtemsq
	res.se2 = se2
	res.se3 = se3
	res.sgh2 = sgh2
	res.sgh3 = sgh3
	res.sgh4 = sgh4
	res.sh2 = sh2
	res.sh3 = sh3
	res.si2 = si2
	res.si3 = si3
	res.sl2 = sl2
	res.sl3 = sl3
	res.sl4 = sl4
	res.s1 = s1
	res.s2 = s2
	res.s3 = s3
	res.s4 = s4
	res.s5 = s5
	res.s6 = s6
	res.s7 = s7
	res.ss1 = ss1
	res.ss2 = ss2
	res.ss3 = ss3
	res.ss4 = ss4
	res.ss5 = ss5
	res.ss6 = ss6
	res.ss7 = ss7
	res.sz1 = sz1
	res.sz2 = sz2
	res.sz3 = sz3
	res.sz11 = sz11
	res.sz12 = sz12
	res.sz13 = sz13
	res.sz21 = sz21
	res.sz22 = sz22
	res.sz23 = sz23
	res.sz31 = sz31
	res.sz32 = sz32
	res.sz33 = sz33
	res.xgh2 = xgh2
	res.xgh3 = xgh3
	res.xgh4 = xgh4
	res.xh2 = xh2
	res.xh3 = xh3
	res.xi2 = xi2
	res.xi3 = xi3
	res.xl2 = xl2
	res.xl3 = xl3
	res.xl4 = xl4
	res.nm = nm
	res.z1 = z1
	res.z2 = z2
	res.z3 = z3
	res.z11 = z11
	res.z12 = z12
	res.z13 = z13
	res.z21 = z21
	res.z22 = z22
	res.z23 = z23
	res.z31 = z31
	res.z32 = z32
	res.z33 = z33
	res.zmol = zmol
	res.zmos = zmos

	return
}
