module Constants

    export me_g, mp_g, me, qe_C, qe, g, eps0

    const me_g = 9.11e-28#9.109383e-28 #g
    const mp_g = 1.*1836.*9.11e-28#1836.152674 * me_g #proton mass in grams
    const me = 1.e-3 * me_g #kg
    const qe_C = 1.6e-19#1.602177e-19
    const qe = 4.8032e-10#4.803204e-10 #statcoulomb
    const g = 5./3. #gamma
    const eps0 = 8.854188e-12

end
