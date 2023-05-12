import pyspedas
from pytplot import tplot


fgm_vars = pyspedas.mms.fgm(trange=['2015-10-16', '2015-10-17'], no_update=no_update)

tplot('mms1_fgm_b_gsm_srvy_l2')
