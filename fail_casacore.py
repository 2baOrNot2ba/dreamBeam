#This code snippet fails with python-casacore version 3.3.1.
#(It works on version 2.2.1)
#Workaround is to call measure a second time.
from casacore.measures import measures
from casacore.quanta import quantity
az=0.1
el=1.0
from_refFrame='J2000'
vr_sph_me = measures().direction(from_refFrame,
                                 quantity(az, 'rad'),
                                 quantity(el, 'rad'))
me = measures()
me.doframe(measures().position('ITRF', '0m', '0m', '0m'))
# Set frame time:
timEpoch = me.epoch('UTC', quantity(12345678, 'd'))
me.do_frame(timEpoch)
to_refFrame='ITRF'
try:
    v_sph_me = me.measure(vr_sph_me, to_refFrame)
except:
    print("measure throw an exception!")
    v_sph_me = me.measure(vr_sph_me, to_refFrame)
print(v_sph_me)

