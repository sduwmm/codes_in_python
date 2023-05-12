from pyspedas import mms_load_fgm
mms_load_fgm(trange=['2015-10-16', '2015-10-17'], data_rate='srvy')
from pyspedas import mms_load_fpi
mms_load_fpi(probe=['3, 4'], trange=['2015-10-16/13:06', '2015-10-16/13:07'],
             data_rate='brst', datatype='des-moms')

from pyspedas.mms import mms_load_mec
mms_load_mec(probe=4)
help(mms_load_mec)
from pyspedas.mms import mms_load_fgm
mms_load_fgm(probe'4', data_rate='brst', trange=['2015-10-16/13:06', '2015-10-16/13:07'],
             time_clip=True)
from pytplot import tplot_names
tplot_names()

from pyspedas.mms import mms_load_edp
mms_load_edp(probe='4', data_rate = 'brst', trange=['2015-10-16/13:06', '2015-10-16/13:07'],
             time_clip=True)
