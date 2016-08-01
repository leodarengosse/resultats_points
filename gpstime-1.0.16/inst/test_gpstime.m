%% Test

[tgps]=ymdhms_t(2010,1,1,12,00,00)

[tgps2]=mjd_t(tgps.mjd)

[tgps3]=jd_t(tgps.jd)

[tgps4]=gpswkd_t(tgps.wk, tgps.wd)

[tgps5]=add_s(tgps,3600)
[tgps6]=add_s(tgps,-86400)

[tgps7]=add_day(tgps,-58)

[tgps8]=gpswks_t(tgps.wk, tgps.wsec)

 [tgps9]=yyyyddds_t(2010,58,56742)

m00(tgps9)

h00(tgps9)

day00(tgps9)

wk00(tgps9)



