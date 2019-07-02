""" Model of integer with nearest neighbor interpolation """

import numpy as np
from myhdl import *

from .utils import IQdata, RTLblocks


class Resampler(object):
    """ Resampler with nearest neighbor interpolation
        List of methods:
        - setRatio - Set ratio and appropriate phase increment;
        - setPhInc - Set phase increment and appropriate ratio;
        - loadFilt - Load filter coefficients;
        - setImpl - Set data format of algorithm implementation;
        - reset - Reset resampler to initial state;
        - loadIntoBuffer - Load samples from d into buffer;
        - step - Do one step in resampling data;
        - rtl - rtl-model of resampler;
        - test_rtl - testbench of rtl-model;
        - convert - Convert myHDL-code to HDL file;
     """

    def __init__(self, ratio, filtCoeffs, phasesNum, dtype=float, impl=('float', None)):
        """
            - ratio - initial resampling ratio, float
            - filtCoeffs - coefficients of filter, np.array
            - phasesNum - number of filter phases
            - dtype - type of input data (and hence output data)
            - impl - format of data for implementation, in [('float', None), ('int', (cSc, phFr)]
                ('float', None) - process data in floating point format
                ('int', (cSc, phFr) - process data in floating point format emulating integer format, where
                    - cSc - width of coefficients in 2'complement, integer
                    - phAccFr - width of phase accumulator for fractional part, integer
                    Filter coefficients in range [-1, 1) is multiplied by 2**(cSc - 1) and rounded to integer format
                    Input data must be integer (or float without fractional part)
        """
        self.iType = None
        self.loadFilt(filtCoeffs, phasesNum)
        self.setImpl(impl)
        self.initRatio = ratio
        self.setRatio(ratio)
        self.dtype = dtype
        self.dataTaps = np.zeros(int(len(self._filt) // self.phasesNum), self.dtype)
        self.dataBuffer = np.array([], dtype=self.dtype)
        self.phase = 0.
        self.wait = False
        self._waitCnt = 0
        self.rtlRefs = np.array([], dtype=self.dtype)
        self.rtlOuts = np.array([], dtype=self.dtype)

    @property
    def ratio(self):
        """ Current resampling ratio """
        return self._ratio

    @property
    def dPh(self):
        """ Current phase increment """
        return self._dPh

    def setRatio(self, val):
        """ Set ratio and appropriate phase increment
            - val - new value of resampling ratio
        """
        self._ratio = float(val)
        if self.iType == 'int':
            self._dPh = np.round((1. / float(val)) * 2 ** self.phFr) / 2 ** self.phFr
        else:
            self._dPh = 1. / float(val)

    def setPhInc(self, val):
        """ Set phase increment and appropriate ratio
            - val - new value of phase increment
        """
        self._ratio = 1. / float(val)
        if self.iType == 'int':
            self._dPh = np.round(float(val) * 2 ** self.phFr) / 2 ** self.phFr
        else:
            self._dPh = float(val)

    def loadFilt(self, filt, phasesNum):
        """ Load filter coefficients.
            - filt - coefficients of lowpass FIR filter, array
            - phasesNum - quantity of phases in polyphase structure, integer
        """
        assert len(filt) % phasesNum == 0, "Number of filter coefficients must be multiple to phases quantity"
        self.phasesNum = phasesNum
        self._filt = filt
        if self.iType is not None:
            self.setImpl((self.iType, self.cSc))

    def setImpl(self, impl):
        """ Set data format of algorithm implementation
            - impl - format of data for implementation, in [('float', None), ('int', (cSc, phFr))]
        """
        assert isinstance(impl, tuple) and len(impl) == 2, (" Impl must be from set [('float', None), "
                                                            "('int', (cSc, phFr))] ")
        self.iType, tmp = impl
        if self.iType == 'int':
            self.cSc, self.phFr = tmp
        else:
            self.cSc, self.phFr = None, None
        if self.iType == 'float':
            self.polyFilt = self._filt.reshape(-1, self.phasesNum).T
        elif (self.iType == 'int') and (isinstance(self.cSc, int)):
            filtInt = np.round(self._filt * 2 ** (self.cSc - 1))
            assert (min(filtInt) >= -2 ** (self.cSc - 1)) and (max(filtInt) < 2 ** (self.cSc - 1)), (
                    " Filter coefficients are " +
                    " too large, they must be in range [-1, (2**(cSc-1) - 0.5) / 2**(cSc-1))] ")
            self.polyFilt = filtInt.reshape(-1, self.phasesNum).T
        else:
            assert False, " Impl must be from set [('float', None), ('int', (cSc, phFr))] "

    def reset(self):
        """ Reset resampler to initial state """
        self.setRatio(self.initRatio)
        self.dataTaps.fill(0.)
        self.dataBuffer = np.array([], dtype=self.dataBuffer.dtype)
        self.phase = 0.
        self.wait = False
        self._waitCnt = 0
        self.rtlRefs = np.array([], dtype=self.dtype)
        self.rtlOuts = np.array([], dtype=self.dtype)

    def _setWait(self, n):
        """ Wait for loading additional n samples in dataBuffer
            n - number of additional required samples
        """
        self._waitCnt = n
        self.wait = True

    def loadIntoBuffer(self, d):
        """ Load samples from d into buffer """
        len_ = 1 if np.isscalar(d) else len(d)
        if len_ >= self._waitCnt:
            self._waitCnt = 0
            self.wait = False
        else:
            self._waitCnt -= len_
        self.dataBuffer = np.append(self.dataBuffer, d)

    def step(self):
        """ Do one step in resampling data"""
        if self.wait:
            return None
        else:
            phInt, phFrac = np.divmod(self.phase, 1)
            phInt = int(phInt)
            if len(self.dataBuffer) < phInt:
                self._setWait(phInt - len(self.dataBuffer))
                return None
            else:
                if phInt > 0:
                    for ii in range(phInt):
                        self.dataTaps = np.roll(self.dataTaps, 1)
                        self.dataTaps[0] = self.dataBuffer[ii]
                    self.dataBuffer = np.delete(self.dataBuffer, np.arange(phInt))
                    self.phase -= phInt
                self.wait = False
                phInd = int(self.phase * self.phasesNum)
                res = np.sum(self.dataTaps * self.polyFilt[phInd, :])
                if self.iType == 'int':
                    if self.dataTaps.dtype == float:
                        res = res // 2 ** (self.cSc - 1)
                    else:  # complex number
                        res = (res.real // 2 ** (self.cSc - 1)) + 1j * (res.imag // 2 ** (self.cSc - 1))
                self.phase += self.dPh
                return res

    @block
    def rtl(self, i_clk, i_rst, i_d, i_dv, i_dPh, i_coeffLd, i_coeff, i_coeffWr,
            i_coeffAdr, o_d, o_dv, o_next):
        """ rtl-model of resampler
            ~~~~ input ports ~~~~
            - i_clk - clock;
            - i_rst - synchronous reset;
            - i_d - input data, for complex input use special class IQdata();
            - i_dv - input data valid;
            - i_dPh - phase increment;
            - i_coeffLd - load coefficients into RAM / resample input signal flag
                True - load coefficients, False - resample;
            - i_coeff - coefficient value;
            - i_coeffWr - write coefficient to RAM (or valid flag for i_coeff and i_coeffAdr ports);
            - i_coeffAdr - address where to write;
            ~~~~ output ports ~~~~
            - o_d - ouput data, see i_d port;
            - o_dv - output data valid;
            - o_next - next data sample request;
        """
        # write to RAM
        phExpand = 1
        intPhLen = phExpand + 2
        nTaps = self.polyFilt.shape[1]
        memWrEn = [Signal(bool(0)) for _ in range(nTaps)]
        coeffIn = Signal(i_coeff.val)
        coeffOut = [Signal(i_coeff.val) for _ in range(nTaps)]
        bank = Signal(intbv(0, min=0, max=nTaps))
        adr = Signal(intbv(0, min=0, max=self.phasesNum))
        assert len(i_coeffAdr) == len(bank) + len(adr), " i_coeffAdr has wrong length "
        adrWr, adrRd = [Signal(adr.val) for _ in range(2)]

        RAMs = [RTLblocks.DualAsRdRAM(i_clk, coeffIn, memWrEn[ii], adrWr, adrRd, coeffOut[ii])
                for ii in range(nTaps)]

        @always_comb
        def adrSlice():
            bank.next = i_coeffAdr[:len(adr)]
            adr.next = i_coeffAdr[len(adr):]

        @always_seq(i_clk.posedge, i_rst)
        def coeffDMX():
            if i_coeffLd:
                for ii in range(len(memWrEn)):
                    if ii == bank:
                        memWrEn[ii].next = i_coeffWr
                    else:
                        memWrEn[ii].next = False
            else:
                for ii in range(len(memWrEn)):
                    memWrEn[ii].next = False
            coeffIn.next = i_coeff
            adrWr.next = adr

        # current phase
        ph, phAcc = [Signal(modbv(0)[len(i_dPh) + phExpand:]) for _ in range(2)]
        nextSample = Signal(bool(0))

        @always_comb
        def phWire():
            ph.next = phAcc + i_dPh

        @always_comb
        def phWire2():
            phMSBs.next = ph[:len(ph) - intPhLen]
            phLSBs.next = ph[len(ph) - intPhLen:]

        phMSBs, cntPh = [Signal(intbv(0)[intPhLen:]) for _ in range(2)]
        phLSBs = Signal(intbv(0)[len(ph) - intPhLen:])
        filtV = Signal(bool(1))
        doutV = Signal(bool(0))

        @always_comb
        def wires():
            adrRd.next = phAcc[len(phAcc) - intPhLen:len(phAcc) - intPhLen - len(adrRd)]
            nextSample.next = phMSBs != 0
            doutV.next = i_dv & filtV

        @always_seq(i_clk.posedge, i_rst)
        def phAccRefresh():
            if i_dv:
                if filtV:
                    if phMSBs > 1:
                        cntPh.next = phMSBs - 2
                        filtV.next = 0
                    else:
                        phAcc.next = concat(intbv(0)[intPhLen:], phLSBs)
                else:
                    if cntPh == 0:
                        filtV.next = 1
                        phAcc.next = concat(intbv(0)[intPhLen:], phLSBs)

        # resample
        if isinstance(i_d, IQdata):
            prodMin = -1 * i_d.I.min * i_coeff.min
            prodMax = i_d.I.max * i_coeff.max
            sumMin = prodMin * nTaps
            sumMax = prodMax * nTaps
            tapsI = [Signal(i_d.I.val) for _ in range(nTaps)]
            tapsQ = [Signal(i_d.I.val) for _ in range(nTaps)]
            prodsI = [Signal(intbv(0, min=prodMin, max=prodMax)) for _ in range(nTaps)]
            prodsQ = [Signal(intbv(0, min=prodMin, max=prodMax)) for _ in range(nTaps)]
            delProdsI = [Signal(prodsI[0].val) for _ in range(nTaps)]
            delProdsQ = [Signal(prodsI[0].val) for _ in range(nTaps)]
            sumsI = [Signal(intbv(0, min=sumMin, max=sumMax)) for _ in range(nTaps - 1)]
            sumsQ = [Signal(intbv(0, min=sumMin, max=sumMax)) for _ in range(nTaps - 1)]

            delays = []

            @always_comb
            def prodBypass():
                delProdsI[0].next = prodsI[0]
                delProdsQ[0].next = prodsQ[0]
                delProdsI[1].next = prodsI[1]
                delProdsQ[1].next = prodsQ[1]

            for ii in range(2, nTaps):
                tmpInI = Signal(prodsI[0].val)
                tmpInI = prodsI[ii]
                tmpInQ = Signal(prodsQ[0].val)
                tmpInQ = prodsQ[ii]
                for jj in range(ii - 1):
                    tmpOutI = Signal(prodsI[0].val)
                    tmpOutQ = Signal(prodsQ[0].val)
                    if jj == ii - 2:
                        delays.append(RTLblocks.Reg(i_clk, i_rst, doutV, tmpInI, delProdsI[ii]))
                        delays.append(RTLblocks.Reg(i_clk, i_rst, doutV, tmpInQ, delProdsQ[ii]))
                    else:
                        delays.append(RTLblocks.Reg(i_clk, i_rst, doutV, tmpInI, tmpOutI))
                        delays.append(RTLblocks.Reg(i_clk, i_rst, doutV, tmpInQ, tmpOutQ))
                        tmpInI = tmpOutI
                        tmpInQ = tmpOutQ

            @always_seq(i_clk.posedge, i_rst)
            def filt():
                if i_dv:
                    for ii in range(nTaps):
                        if nextSample:
                            if ii == 0:
                                tapsI[ii].next = i_d.I
                                tapsQ[ii].next = i_d.Q
                            else:
                                tapsI[ii].next = tapsI[ii - 1]
                                tapsQ[ii].next = tapsQ[ii - 1]
                if doutV:
                    for ii in range(nTaps):
                        prodsI[ii].next = coeffOut[ii] * tapsI[ii]
                        prodsQ[ii].next = coeffOut[ii] * tapsQ[ii]
                    for ii in range(nTaps - 1):
                        if ii == 0:
                            sumsI[ii].next = delProdsI[ii] + delProdsI[ii + 1]
                            sumsQ[ii].next = delProdsQ[ii] + delProdsQ[ii + 1]
                        else:
                            sumsI[ii].next = sumsI[ii - 1] + delProdsI[ii + 1]
                            sumsQ[ii].next = sumsQ[ii - 1] + delProdsQ[ii + 1]

            @always_comb
            def outputD():
                o_d.I.next = sumsI[nTaps - 2] >> (len(i_coeff) - 1)
                o_d.Q.next = sumsQ[nTaps - 2] >> (len(i_coeff) - 1)
        else:
            prodMin = -1 * i_d.min * i_coeff.min
            prodMax = i_d.max * i_coeff.max
            sumMin = prodMin * nTaps
            sumMax = prodMax * nTaps
            taps = [Signal(i_d.val) for _ in range(nTaps)]
            prods = [Signal(intbv(0, min=prodMin, max=prodMax)) for _ in range(nTaps)]
            delProds = [Signal(prods[0].val) for _ in range(nTaps)]
            sums = [Signal(intbv(0, min=sumMin, max=sumMax)) for _ in range(nTaps - 1)]

            delays = []

            @always_comb
            def prodBypass():
                delProds[0].next = prods[0]
                delProds[1].next = prods[1]

            for ii in range(2, nTaps):
                tmpIn = Signal(prods[0].val)
                tmpIn = prods[ii]
                for jj in range(ii - 1):
                    tmpOut = Signal(prods[0].val)
                    if jj == ii - 2:
                        delays.append(RTLblocks.Reg(i_clk, i_rst, doutV, tmpIn, delProds[ii]))
                    else:
                        delays.append(RTLblocks.Reg(i_clk, i_rst, doutV, tmpIn, tmpOut))
                        tmpIn = tmpOut

            @always_seq(i_clk.posedge, i_rst)
            def filt():
                if i_dv:
                    for ii in range(nTaps):
                        if nextSample:
                            if ii == 0:
                                taps[ii].next = i_d
                            else:
                                taps[ii].next = taps[ii - 1]
                if doutV:
                    for ii in range(nTaps):
                        prods[ii].next = coeffOut[ii] * taps[ii]
                    for ii in range(nTaps - 1):
                        if ii == 0:
                            sums[ii].next = delProds[ii] + delProds[ii + 1]
                        else:
                            sums[ii].next = sums[ii - 1] + delProds[ii + 1]

            @always_comb
            def outputD():
                o_d.next = sums[nTaps - 2] >> (len(i_coeff) - 1)

        # output
        @always_comb
        def output():
            o_next.next = nextSample
            o_dv.next = doutV

        return instances()

    @block
    def test_rtl(self, d, dW, period=10, wrToFile=False):
        """ Testbench. Compare output of rtl-model with 'int' implementation
            - d - stimuli for models, array of integer;
            - dW - input/output data bitwidth;
            - period of clock, integer;
            - wrToFile - if True, then stimuli and output write to file
        """
        assert self.iType == 'int', (" Testbench of rtl-model is valid only for 'int' implementation. " +
                                     " But current implementation is %s " % self.iType)
        assert self.ratio > 0.4, (" Rtl-model works only for ratio greater then 2/5, current ratio is " +
                                  f"{self.ratio}")
        tapsN = self.polyFilt.shape[1]
        delayRTL = tapsN - 1
        bankLen = int(np.ceil(np.log2(tapsN)))
        RAMLen = int(np.ceil(np.log2(self.polyFilt.shape[0])))
        adrLen = RAMLen + bankLen
        i_clk = Signal(bool(0))
        i_rst = ResetSignal(0, active=bool(1), isasync=False)
        if self.dtype == float:
            i_d = Signal(intbv(0, min=-2 ** (dW - 1), max=2 ** (dW - 1)))
            o_d = Signal(intbv(0, min=-2 ** (dW - 1), max=2 ** (dW - 1)))
        else:
            i_d = IQdata(dW)
            o_d = IQdata(dW)
        i_dv = Signal(bool(0))
        i_dPh = Signal(intbv(0)[32:])

        i_coeffLd = Signal(bool(0))
        i_coeff = Signal(intbv(0, min=-2 ** (self.cSc - 1), max=2 ** (self.cSc - 1)))
        i_coeffWr = Signal(bool(0))
        i_coeffAdr = Signal(intbv(0, min=0, max=2 ** adrLen))
        o_dv = Signal(bool(0))
        o_next = Signal(bool(0))

        uut = self.rtl(i_clk, i_rst, i_d, i_dv, i_dPh, i_coeffLd, i_coeff, i_coeffWr,
                       i_coeffAdr, o_d, o_dv, o_next)

        @always(delay(int(period // 2)))
        def clk_gen():
            i_clk.next = not i_clk

        if wrToFile:
            inFile = open('./vhdl/stim/iStimuli_Resampler.txt', 'w')

            if self.dtype == float:
                @instance
                def wrInStimuli():
                    while True:
                        if self.dtype == float:
                            yield i_clk.posedge
                            print((f"{int(i_rst.val)} {int(i_coeffLd.val)} {int(i_coeffWr.val)} {i_coeffAdr.val} " +
                                   f"{i_coeff.val} {i_d.val} {int(i_dv.val)} {i_dPh.val}"), file=inFile)
                            yield i_clk.negedge
                            print((f"{int(i_rst.val)} {int(i_coeffLd.val)} {int(i_coeffWr.val)} {i_coeffAdr.val} " +
                                   f"{i_coeff.val} {i_d.val} {int(i_dv.val)} {i_dPh.val}"), file=inFile)
            else:
                @instance
                def wrInStimuli():
                    while True:
                        if self.dtype == float:
                            yield i_clk.posedge
                        yield i_clk.posedge
                        print((f"{int(i_rst.val)} {int(i_coeffLd.val)} {int(i_coeffWr.val)} {i_coeffAdr.val} " +
                               f"{i_coeff.val} {i_d.I.val} {i_d.Q.val} {int(i_dv.val)} {i_dPh.val}"), file=inFile)
                        yield i_clk.negedge
                        print((f"{int(i_rst.val)} {int(i_coeffLd.val)} {int(i_coeffWr.val)} {i_coeffAdr.val} " +
                               f"{i_coeff.val} {i_d.I.val} {i_d.Q.val} {int(i_dv.val)} {i_dPh.val}"), file=inFile)

            outFile = open('./vhdl/stim/oStimuli_Resampler.txt', 'w')

            if self.dtype == float:
                @always(i_clk.posedge)
                def wrOut():
                    if o_dv:
                        print(f"{o_d.val} {int(o_next.val)}", file=outFile)
            else:
                @always(i_clk.posedge)
                def wrOut():
                    if o_dv:
                        print(f"{o_d.I.val} {o_d.Q.val} {int(o_next.val)}", file=outFile)

        if self.dtype == float:
            @instance
            def stimuli():
                yield i_clk.posedge
                i_rst.next = bool(1)
                for _ in range(3):
                    yield i_clk.posedge
                i_rst.next = bool(0)

                # load coefficients to filter
                yield i_clk.negedge
                i_coeffLd.next = bool(1)
                i_coeffWr.next = bool(1)
                nPh, nTaps = self.polyFilt.shape
                for ii in range(nTaps):
                    for jj in range(nPh):
                        i_coeffAdr.next[:RAMLen] = ii
                        i_coeffAdr.next[RAMLen:] = jj
                        i_coeff.next = int(self.polyFilt[jj, ii])
                        yield i_clk.negedge
                i_coeffWr.next = bool(0)
                yield i_clk.negedge
                i_coeffLd.next = bool(0)
                i_rst.next = bool(1)
                yield i_clk.negedge
                i_rst.next = bool(0)

                for _ in range(3):
                    yield i_clk.negedge
                nTick = 0
                self.reset()

                self.loadIntoBuffer(d)
                val = 0
                while val is not None:
                    val = self.step()
                    if val is not None:
                        self.rtlRefs = np.append(self.rtlRefs, val)

                nData = 0
                refInd = 0
                i_dPh.next = int(np.round(self.dPh * 2 ** 30))
                if self.ratio > 1.:
                    drtl = np.append(0., d)
                else:
                    drtl = d
                while nData < len(drtl):
                    i_d.next = int(drtl[nData])
                    i_dv.next = bool(1)
                    yield i_clk.negedge
                    if o_next:
                        nData += 1
                    if o_dv:
                        rtlOut = int(o_d.val)
                        self.rtlOuts = np.append(self.rtlOuts, rtlOut)
                        if nTick < delayRTL:
                            nTick += 1
                        else:
                            assert rtlOut == self.rtlRefs[refInd], (f" Outputs from rtl-model and reference " +
                                                                    f" doesn't match, {rtlOut} != " +
                                                                    f"{self.rtlRefs[refInd]}, refInd is {refInd} ")
                            refInd += 1
                raise StopSimulation
        else:
            @instance
            def stimuli():
                yield i_clk.posedge
                i_rst.next = bool(1)
                for _ in range(3):
                    yield i_clk.posedge
                i_rst.next = bool(0)

                # load coefficients to filter
                yield i_clk.negedge
                i_coeffLd.next = bool(1)
                i_coeffWr.next = bool(1)
                nPh, nTaps = self.polyFilt.shape
                for ii in range(nTaps):
                    for jj in range(nPh):
                        i_coeffAdr.next[:RAMLen] = ii
                        i_coeffAdr.next[RAMLen:] = jj
                        i_coeff.next = int(self.polyFilt[jj, ii])
                        yield i_clk.negedge
                i_coeffWr.next = bool(0)
                yield i_clk.negedge
                i_coeffLd.next = bool(0)
                i_rst.next = bool(1)
                yield i_clk.negedge
                i_rst.next = bool(0)

                for _ in range(3):
                    yield i_clk.negedge

                nTick = 0
                self.reset()

                self.loadIntoBuffer(d)
                val = 0
                while val is not None:
                    val = self.step()
                    if val is not None:
                        self.rtlRefs = np.append(self.rtlRefs, val)

                nData = 0
                refInd = 0
                i_dPh.next = int(np.round(self.dPh * 2 ** 30))
                if self.ratio > 1.:
                    drtl = np.append(0. + 1j * 0., d)
                else:
                    drtl = d
                while nData < len(drtl):
                    i_d.I.next = int(drtl[nData].real)
                    i_d.Q.next = int(drtl[nData].imag)
                    i_dv.next = bool(1)
                    yield i_clk.negedge
                    if o_next:
                        nData += 1
                    if o_dv:
                        rtlOut = int(o_d.I.val) + 1j * int(o_d.Q.val)
                        self.rtlOuts = np.append(self.rtlOuts, rtlOut)
                        if nTick < delayRTL:
                            nTick += 1
                        else:
                            assert rtlOut == self.rtlRefs[refInd], (f" Outputs from rtl-model and reference " +
                                                                    f" doesn't match, {rtlOut} != " +
                                                                    f"{self.rtlRefs[refInd]}, refInd is {refInd} ")
                            refInd += 1
                raise StopSimulation

        return instances()

    def convert(self, dW, hdl='VHDL', name='Resampler', path='./vhdl'):
        """ Convert myHDL-code to HDL file.
            - dW - input/output data bitwidth
            - hdl - HDL language in ['VHDL', 'Verilog']
            - name - name of the entity (module) and appropriate file, string
            - path - destination folder, string
        """
        bankLen = int(np.ceil(np.log2(self.polyFilt.shape[1])))
        adrLen = int(np.ceil(np.log2(self.polyFilt.shape[0])) + bankLen)
        i_clk = Signal(bool(0))
        i_rst = ResetSignal(0, active=bool(1), isasync=False)
        if self.dtype == float:
            i_d = Signal(intbv(0, min=-2 ** (dW - 1), max=2 ** (dW - 1)))
            o_d = Signal(intbv(0, min=-2 ** (dW - 1), max=2 ** (dW - 1)))
        else:
            i_d = IQdata(dW)
            o_d = IQdata(dW)
        i_dv = Signal(bool(0))
        i_dPh = Signal(intbv(0)[32:])

        i_coeffLd = Signal(bool(0))
        i_coeff = Signal(intbv(0, min=-2 ** (self.cSc - 1), max=2 ** (self.cSc - 1)))
        i_coeffWr = Signal(bool(0))
        i_coeffAdr = Signal(intbv(0, min=0, max=2 ** adrLen))
        o_dv = Signal(bool(0))
        o_next = Signal(bool(0))

        inst = self.rtl(i_clk, i_rst, i_d, i_dv, i_dPh, i_coeffLd, i_coeff, i_coeffWr,
                        i_coeffAdr, o_d, o_dv, o_next)

        inst.convert(hdl=hdl, name=name, path=path)