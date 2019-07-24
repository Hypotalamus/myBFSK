""" Run all test using this class """

import numpy as np
import scipy.signal as signal
import pytest

from src.Decimator import Decimator
from src.NCO import NCO
from src.Resampler import Resampler
from src.CORDIC import CORDIC

class TestThemAll(object):
    """ Collection of tests """

    def test_Decimator(self, sigType=complex):
        """ Run test of polyphase decimator. Some decimator with factor 4 is used.
            - sigType - data type of input data in [float, complex]
        """

        decF = 4
        nFilt = decF * 16
        filt = signal.firwin2(nFilt, [0., 1/4/decF, 1/2/decF, 1.], [1., 1., 0., 0.])
        filt = filt / np.sum(filt)
        decImpl = ('int', 16)

        dec = Decimator(filt, decF, dtype=sigType, impl=decImpl)

        if sigType == float:
            sig = np.cos(2 * np.pi * .05 / 16 * np.arange(5000)) + .5 * np.cos(
                2 * np.pi * 3 * .05 / 16 * np.arange(5000))
        else:
            sig = np.exp(1j * 2 * np.pi * .05 / 16 * np.arange(5000)) + .5 * np.exp(
                1j * 2 * np.pi * 3 * .05 / 16 * np.arange(5000))

        sigSc = 12
        sigNorm =  .99 * sig / max(sig) * (2**(sigSc-1)-1) / (2**(sigSc-1))
        sigInt = np.round(sigNorm * 2**(sigSc-1))

        decFmax = 16

        tb = dec.test_rtl(sigInt, sigSc, decFmax)
        tb.config_sim(trace=False)
        tb.run_sim()

    def test_NCO(self, f0=216, fs=10000, ampl=np.sqrt(2.), phi0=.1,
                 impl=('int', (32, 10, 12))):
        """ Run test of quadrature Numeric Controlled Oscillator.
        :param f0: frequency of generated complex exponential;
        :param fs: sampling frequency;
        :param ampl: amplitude of cosine/sine waveforms;
        :param phi0: initial phase in [0, 1];
        :param impl: implementation, 'int' (integer) is required, 2nd
            argument is (phW, ROMadrW, outW), where
            - phW - bitwidth of NCO phase, integer > 0
            - ROMadrW - bitwidth of ROM address with sine/cosine samples (1/4-period ROM)
            - outW - bitwidth of NCO output sample, integer > 0
        :return:
        """
        nco = NCO(f0, fs, ampl=ampl, phi0=phi0, impl=impl)
        tb = nco.test_rtl()
        tb.config_sim(trace=False)
        tb.run_sim()

    @pytest.mark.parametrize("ratio", {5/9, 1., np.sqrt(2)})
    def test_Resampler(self, ratio, sigType=complex):
        """ Run test of resampler.
            - ratio - resampling ratio;
            - sigType - data type of input data in [float, complex]
        """

        N = 4096
        taps = 16
        nFilt = N * taps
        filt = signal.firwin2(nFilt, [0., .35 / N, .5 / N, 1.], [1., 1., 0., 0.])
        filt = N * filt / np.sum(filt)
        impl = ('int', (16, 30))

        res = Resampler(ratio, filt, N, dtype=sigType, impl=impl)

        if sigType == float:
            sig = (np.cos(2 * np.pi * 1/50 * np.arange(1000)) + 1.2 * np.cos(
                2 * np.pi * 3 * 1/50 * np.arange(1000)) -
                   .3 * np.cos( 2 * np.pi * 3.5 * 1/50 * np.arange(1000)))
        else:
            sig = (np.exp(1j * 2 * np.pi * 1/50 * np.arange(1000)) + 1.2 * np.exp(
                1j * 2 * np.pi * 3 * 1/50 * np.arange(1000)) -
                   .3 * np.exp(1j * 2 * np.pi * 3.5 * 1/50 * np.arange(1000)))

        sigSc = 12
        sigNorm =  .99 * sig / max(sig) * (2**(sigSc-1)-1) / (2**(sigSc-1))
        sigInt = np.round(sigNorm * 2**(sigSc-1))

        tb = res.test_rtl(sigInt, sigSc)
        tb.config_sim(trace=False)
        tb.run_sim()

    def test_CORDIC(self, dBw=12, cordBw=16, cordBwInt=19, stQnt=17,
                    rndQnt=int(1e4)):
        """ Run test of CORDIC processor
            - dBw - input data (I and Q) bitwidth;
            - cordBw - output angle bitwidth;
            - cordBwInt - internal data bitwidth;
            - stQnt - number of rotational stages;
            - rndQnt - number of test random data
        """
        cord = CORDIC(dBw, cordBw, stQnt, cordBwInt)
        tb = cord.test_rtl(randQnt=rndQnt)
        tb.config_sim(trace=False)
        tb.run_sim()