""" Model of integer decimator via polyphase FIR filter """

import numpy as np
from myhdl import *

from .utils import IQdata, RTLblocks

class Decimator(object):
    """ Model of integer decimator via polyphase FIR filter
        List of methods:
        - loadFilt - Load filter coefficients;
        - setImpl - Set data format of algorithm implementation ('float' or 'fixed-point')
        - reset - Reset decimator;
        - decimate - Decimate input signal step-by-step;
        - rtl - RTL-model of polyphase FIR filter in myHDL;
        - test_rtl - Testbench. Compare output of rtl-model with 'int' implementation;
        - convert - Convert myHDL-code to HDL file;
    """

    def __init__(self, filt, decF, dtype=float, impl=('float', None)):
        """
            - filt - coefficients of lowpass FIR filter, array
            - decF - decimation factor (number of phases); len(filt)/decF must be integer, integer
            - dtype - type of input data (and hence output data)
            - impl - format of data for implementation, in [('float', None), ('int', cSc)]
                ('float', None) - process data in floating point format
                ('int', cSc) - process data in floating point format emulating integer format, where
                    - cSc - width of coefficients in 2'complement, integer
                    Filter coefficients in range [-1, 1) is multiplied by 2**(cSc - 1) and rounded to integer format
                    Input data must be integer (or float without fractional part)
        """
        self.iType = None
        self.loadFilt(filt, decF)
        self.setImpl(impl)
        self.dtype = dtype
        self.dataTaps = np.zeros(self.polyFilt.shape, self.dtype)
        self.phSel = 0
        self.rtlRefs = np.array([], dtype=self.dtype)
        self.rtlOuts = np.array([], dtype=self.dtype)

    def loadFilt(self, filt, decF):
        """ Load filter coefficients.
            - filt - coefficients of lowpass FIR filter, array
            - decF - decimation factor (number of phases); len(filt)/decF must be integer, integer
        """
        assert len(filt) % decF == 0, " Length of filter must be multiple of decimation factor "
        self.decF = decF
        self._filt = filt
        if self.iType is not None:
            self.setImpl((self.iType, self.cSc))

    def setImpl(self, impl):
        """ Set data format of algorithm implementation
            - impl - format of data for implementation, in [('float', None), ('int', cSc)]
        """
        assert isinstance(impl, tuple) and len(impl) == 2, " Impl must be from set [('float', None), ('int', cSc)] "
        self.iType, self.cSc = impl
        if self.iType == 'float':
            self.polyFilt = self._filt.reshape(-1, self.decF).T
        elif (self.iType == 'int') and (isinstance(self.cSc, int)):
            filtInt = np.round(self._filt * 2 ** (self.cSc - 1))
            assert (min(filtInt) >= -2 ** (self.cSc - 1)) and (max(filtInt) < 2 ** (self.cSc - 1)), (
                        " Filter coefficients are " +
                        " too large, they must be in range [-1, (2**(cSc-1) - 0.5) / 2**(cSc-1))] ")
            self.polyFilt = filtInt.reshape(-1, self.decF).T
        else:
            assert False, " Impl must be from set [('float', None), ('int', cSc)] "

    def reset(self):
        """ Reset decimator """
        self.dataTaps.fill(0.)
        self.phSel = 0
        self.rtlRefs = np.array([], dtype=self.dtype)
        self.rtlOuts = np.array([], dtype=self.dtype)

    def decimate(self, d):
        """ Decimate input signal
            - d - input sample; if self.iType = 'int' then input data must be integer
                    (or float without fractional part)
            Return output sample every decFth input sample, None otherwise
        """
        if self.iType == 'int':
            assert (d.real % 1. == 0.) and (d.imag % 1. == 0.), ("Error. Input data must be " +
                                                                 " integer (or float without fractional part) ")
        self.dataTaps[self.phSel, :-1] = self.dataTaps[self.phSel, 1:]
        self.dataTaps[self.phSel, -1] = d
        if self.phSel == self.decF - 1:
            self.phSel = 0
            res = np.sum(self.dataTaps * self.polyFilt)
            if self.iType == 'int':
                if self.dataTaps.dtype == float:
                    res = res // 2 ** (self.cSc - 1)
                else:  # complex number
                    res = (res.real // 2 ** (self.cSc - 1)) + 1j * (res.imag // 2 ** (self.cSc - 1))
            return res
        else:
            self.phSel += 1
            return None

    @block
    def rtl(self, i_clk, i_rst, i_d, i_dv, i_decF, i_coeffLd, i_coeff, i_coeffWr,
            i_coeffAdr, o_d, o_dv, nTaps, decFmax):
        """ RTL-model of polyphase FIR filter in myHDL
            ~~~~ input ports ~~~~
            - i_clk - clock;
            - i_rst - synchronous reset (RAM with coefficients doesn't clear though);
            - i_d - input data, for complex input use special class IQdata();
            - i_dv - input data valid flag;
            - i_decF - current decimation factor;
            - i_coeffLd - load coefficients into RAM / decimate input signal flag
                True - load coefficients, False - decimate;
            - i_coeff - coefficient value;
            - i_coeffWr - write coefficient to RAM (or valid flag for i_coeff and i_coeffAdr ports);
            - i_coeffAdr - address where to write;
            ~~~~ output ports ~~~~
            - o_d - output data, see i_d port;
            - o_dv - output data valid flag;
            ~~~~ parameters ~~~~
            - nTaps - quantity of taps in filter (quantity of phases). Overall length of filter must be nTaps * i_decF;
            - decFmax - maximum value of decimation factor
        """
        # write to RAM
        memWrEn = [Signal(bool(0)) for _ in range(nTaps)]
        coeffIn = Signal(i_coeff.val)
        coeffOut = [Signal(i_coeff.val) for _ in range(nTaps)]
        cBank = Signal(intbv(0, min=0, max=nTaps))
        cAdr = Signal(intbv(0, min=0, max=decFmax))
        assert len(i_coeffAdr) == len(cBank) + len(cAdr), " i_coeffAdr has wrong length "
        cAdrReg = Signal(cAdr.val)
        cAdrReg2 = Signal(cAdr.val)

        RAMs = [RTLblocks.SingleAsRdRAM(i_clk, coeffIn, memWrEn[ii], cAdrReg, coeffOut[ii])
                for ii in range(nTaps)]

        @always_comb
        def adrSlice():
            cBank.next = i_coeffAdr[:len(cAdr)]
            cAdr.next = i_coeffAdr[len(cAdr):]

        @always(i_clk.posedge)
        def coeffDMX():
            if i_rst:
                for ii in range(len(memWrEn)):
                    memWrEn[ii].next = False
                coeffIn.next = 0
                cAdrReg.next = 0
                cAdrReg2.next = 0
            else:
                if i_coeffLd:
                    for ii in range(len(memWrEn)):
                        if ii == cBank:
                            memWrEn[ii].next = i_coeffWr
                        else:
                            memWrEn[ii].next = False
                    coeffIn.next = i_coeff
                    cAdrReg.next = cAdr
                else:
                    for ii in range(len(memWrEn)):
                        memWrEn[ii].next = False
                    coeffIn.next = 0
                    if i_dv:
                        if cAdrReg >= i_decF:
                            cAdrReg.next = 0
                        else:
                            cAdrReg.next = cAdrReg + 1
                cAdrReg2.next = cAdrReg

        # decimate input data

        accClr = Signal(bool(0))
        filtV = Signal(bool(0))

        if isinstance(i_d, IQdata):
            dI, dQ = [Signal(intbv(i_d.I.val)) for _ in range(2)]

            @always_comb
            def inData():
                dI.next = i_d.I
                dQ.next = i_d.Q

            MACoutMin = -1 * i_coeff.min * i_d.I.min
            MACoutMax = i_coeff.max * i_d.I.max
            fGain = int(np.ceil(np.sum(np.abs(self._filt))))
            wAcc = int(np.ceil(np.log2(fGain * MACoutMax)))
            MACoutI = [Signal(intbv(0, min=fGain * MACoutMin, max=fGain * MACoutMax)) for _ in range(nTaps)]
            MACoutQ = [Signal(intbv(0, min=fGain * MACoutMin, max=fGain * MACoutMax)) for _ in range(nTaps)]

            MACsI = [RTLblocks.MAC(i_clk, i_rst, dI, coeffOut[ii], i_dv, accClr, MACoutI[ii], None, None,
                                   wAcc=wAcc)
                     for ii in range(nTaps)]
            MACsQ = [RTLblocks.MAC(i_clk, i_rst, dQ, coeffOut[ii], i_dv, accClr, MACoutQ[ii], None, None,
                                   wAcc=wAcc)
                     for ii in range(nTaps)]

            filtRegsI = [Signal(intbv(0, min=fGain * MACoutMin, max=fGain * MACoutMax)) for _ in range(nTaps)]
            filtRegsQ = [Signal(intbv(0, min=fGain * MACoutMin, max=fGain * MACoutMax)) for _ in range(nTaps)]

            @always(i_clk.posedge)
            def filt():
                if i_rst:
                    for ii in range(nTaps):
                        filtRegsI[ii].next = 0
                        filtRegsQ[ii].next = 0
                    accClr.next = False
                    filtV.next = False
                else:
                    if i_dv and (cAdrReg == i_decF):
                        accClr.next = True
                    else:
                        accClr.next = False
                    if i_dv and (cAdrReg2 == i_decF):
                        for ii in range(nTaps):
                            if ii == 0:
                                filtRegsI[ii].next = MACoutI[ii]
                                filtRegsQ[ii].next = MACoutQ[ii]
                            else:
                                filtRegsI[ii].next = MACoutI[ii] + filtRegsI[ii - 1]
                                filtRegsQ[ii].next = MACoutQ[ii] + filtRegsQ[ii - 1]
                        filtV.next = True
                    else:
                        filtV.next = False

            oMax = len(i_d.I) + self.cSc - 1
            oMin = self.cSc - 1

            @always_comb
            def output():
                o_d.I.next = filtRegsI[nTaps - 1][oMax:oMin].signed()
                o_d.Q.next = filtRegsQ[nTaps - 1][oMax:oMin].signed()
                o_dv.next = filtV and i_dv

        else:
            MACoutMin = -1 * i_coeff.min * i_d.min
            MACoutMax = i_coeff.max * i_d.max
            fGain = int(np.ceil(np.sum(np.abs(self._filt))))
            wAcc = int(np.ceil(np.log2(fGain * MACoutMax)))
            MACout = [Signal(intbv(0, min=fGain * MACoutMin, max=fGain * MACoutMax)) for _ in range(nTaps)]

            MACs = [RTLblocks.MAC(i_clk, i_rst, i_d, coeffOut[ii], i_dv, accClr, MACout[ii], None, None,
                                  wAcc=wAcc)
                    for ii in range(nTaps)]

            filtRegs = [Signal(intbv(0, min=fGain * MACoutMin, max=fGain * MACoutMax)) for _ in range(nTaps)]

            @always(i_clk.posedge)
            def filt():
                if i_rst:
                    for ii in range(nTaps):
                        filtRegs[ii].next = 0
                    accClr.next = False
                    filtV.next = False
                else:
                    if i_dv and (cAdrReg == i_decF):
                        accClr.next = True
                    else:
                        accClr.next = False
                    if i_dv and (cAdrReg2 == i_decF):
                        for ii in range(nTaps):
                            if ii == 0:
                                filtRegs[ii].next = MACout[ii]
                            else:
                                filtRegs[ii].next = MACout[ii] + filtRegs[ii - 1]
                        filtV.next = True
                    else:
                        filtV.next = False

            oMax = len(i_d) + self.cSc - 1
            oMin = self.cSc - 1

            @always_comb
            def output():
                o_d.next = filtRegs[nTaps - 1][oMax:oMin].signed()
                o_dv.next = filtV and i_dv

        return instances()

    @block
    def test_rtl(self, d, dW, decFmax, period=10, wrToFile=False):
        """ Testbench. Compare output of rtl-model with 'int' implementation
            - d - stimuli for models, array of integer;
            - dW - input/output data bitwidth;
            - decFmax - maximum decimation ratio for rtl;
            - period of clock, integer;
            - wrToFile - if True, then stimuli and output write to file
        """
        assert self.iType == 'int', (" Testbench of rtl-model is valid only for 'int' implementation. " +
                                     " But current implementation is %s " % self.iType)
        delayRTL = 1
        nTaps = self.polyFilt.shape[1]
        cBankLen = int(np.ceil(np.log2(decFmax)))
        cAdrLen = int(np.ceil(np.log2(nTaps)) + cBankLen)
        i_clk = Signal(bool(0))
        i_rst = Signal(bool(0))
        if self.dtype == float:
            i_d = Signal(intbv(0, min=-2 ** (dW - 1), max=2 ** (dW - 1)))
            o_d = Signal(intbv(0, min=-2 ** (dW - 1), max=2 ** (dW - 1)))
        else:
            i_d = IQdata(dW)
            o_d = IQdata(dW)
        i_dv = Signal(bool(0))
        i_decF = Signal(intbv(self.decF - 1, min=0, max=decFmax))
        i_coeffLd = Signal(bool(0))
        i_coeff = Signal(intbv(0, min=-2 ** (self.cSc - 1), max=2 ** (self.cSc - 1)))
        i_coeffWr = Signal(bool(0))
        i_coeffAdr = Signal(intbv(0, min=0, max=2 ** cAdrLen))
        o_dv = Signal(bool(0))

        uut = self.rtl(i_clk, i_rst, i_d, i_dv, i_decF, i_coeffLd, i_coeff, i_coeffWr,
                       i_coeffAdr, o_d, o_dv, nTaps, decFmax)

        @always(delay(int(period // 2)))
        def clk_gen():
            i_clk.next = not i_clk

        if wrToFile:
            inFile = open('./vhdl/stim/iStimuli_Decimator.txt', 'w')

            if self.dtype == float:
                @instance
                def wrInStimuli():
                    while True:
                        if self.dtype == float:
                            yield i_clk.posedge
                            print((f"{int(i_rst.val)} {int(i_coeffLd.val)} {int(i_coeffWr.val)} {i_coeffAdr.val} " +
                                   f"{i_coeff.val} {i_d.val} {int(i_dv.val)} {i_decF.val}"), file=inFile)
                            yield i_clk.negedge
                            print((f"{int(i_rst.val)} {int(i_coeffLd.val)} {int(i_coeffWr.val)} {i_coeffAdr.val} " +
                                   f"{i_coeff.val} {i_d.val} {int(i_dv.val)} {i_decF.val}"), file=inFile)
            else:
                @instance
                def wrInStimuli():
                    while True:
                        if self.dtype == float:
                            yield i_clk.posedge
                        yield i_clk.posedge
                        print((f"{int(i_rst.val)} {int(i_coeffLd.val)} {int(i_coeffWr.val)} {i_coeffAdr.val} " +
                               f"{i_coeff.val} {i_d.I.val} {i_d.Q.val} {int(i_dv.val)} {i_decF.val}"), file=inFile)
                        yield i_clk.negedge
                        print((f"{int(i_rst.val)} {int(i_coeffLd.val)} {int(i_coeffWr.val)} {i_coeffAdr.val} " +
                               f"{i_coeff.val} {i_d.I.val} {i_d.Q.val} {int(i_dv.val)} {i_decF.val}"), file=inFile)

            outFile = open('./vhdl/stim/oStimuli_Decimator.txt', 'w')

            if self.dtype == float:
                @always(i_clk.posedge)
                def wrOut():
                    if o_dv:
                        print(f"{o_d.val}", file=outFile)
            else:
                @always(i_clk.posedge)
                def wrOut():
                    if o_dv:
                        print(f"{o_d.I.val} {o_d.Q.val}", file=outFile)

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
                        i_coeffAdr.next[:cBankLen] = nTaps - 1 - ii
                        i_coeffAdr.next[cBankLen:] = nPh - 1 - jj
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
                self.rtlRefs = np.array([], dtype=self.dtype)
                self.rtlOuts = np.array([], dtype=self.dtype)
                for din in d:
                    ref = self.decimate(din)
                    if ref:
                        self.rtlRefs = np.append(self.rtlRefs, ref)
                    i_d.next = int(din)
                    i_dv.next = bool(1)
                    yield i_clk.negedge
                    if o_dv:
                        rtlOut = int(o_d.val)
                        self.rtlOuts = np.append(self.rtlOuts, rtlOut)
                        if nTick < delayRTL:
                            nTick += 1
                        else:
                            assert rtlOut == self.rtlRefs[-1], (" Outputs from rtl-model and reference " +
                                                                "doesn't match ")
                print(' **** Test passed **** ')
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
                        i_coeffAdr.next[:cBankLen] = nTaps - 1 - ii
                        i_coeffAdr.next[cBankLen:] = nPh - 1 - jj
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
                self.rtlRefs = np.array([], dtype=self.dtype)
                self.rtlOuts = np.array([], dtype=self.dtype)
                for din in d:
                    ref = self.decimate(din)
                    if ref:
                        self.rtlRefs = np.append(self.rtlRefs, ref)
                    i_d.I.next = int(din.real)
                    i_d.Q.next = int(din.imag)
                    i_dv.next = bool(1)
                    yield i_clk.negedge
                    if o_dv:
                        rtlOut = int(o_d.I.val) + 1j * int(o_d.Q.val)
                        self.rtlOuts = np.append(self.rtlOuts, rtlOut)
                        if nTick < delayRTL:
                            nTick += 1
                        else:
                            assert rtlOut == self.rtlRefs[-1], (" Outputs from rtl-model and reference " +
                                                                "doesn't match ")
                print(' **** Test passed **** ')
                raise StopSimulation

        return instances()

    def convert(self, dW, decFmax, hdl='VHDL', name='PolyphaseDecimator', path='./vhdl'):
        """ Convert myHDL-code to HDL file.
            - dW - input/output data bitwidth
            - decFmax - maximum decimation ratio
            - hdl - HDL language in ['VHDL', 'Verilog']
            - name - name of the entity (module) and appropriate file, string
            - path - destination folder, string
        """
        nTaps = self.polyFilt.shape[1]
        cBankLen = int(np.ceil(np.log2(decFmax)))
        cAdrLen = int(np.ceil(np.log2(nTaps)) + cBankLen)
        i_clk = Signal(bool(0))
        i_rst = Signal(bool(0))
        if self.dtype == float:
            i_d = Signal(intbv(0, min=-2 ** (dW - 1), max=2 ** (dW - 1)))
            o_d = Signal(intbv(0, min=-2 ** (dW - 1), max=2 ** (dW - 1)))
        else:
            i_d = IQdata(dW)
            o_d = IQdata(dW)
        i_dv = Signal(bool(0))
        i_decF = Signal(intbv(0, min=0, max=decFmax))
        i_coeffLd = Signal(bool(0))
        i_coeff = Signal(intbv(0, min=-2 ** (self.cSc - 1), max=2 ** (self.cSc - 1)))
        i_coeffWr = Signal(bool(0))
        i_coeffAdr = Signal(intbv(0, min=0, max=2 ** cAdrLen))
        o_dv = Signal(bool(0))

        inst = self.rtl(i_clk, i_rst, i_d, i_dv, i_decF, i_coeffLd, i_coeff, i_coeffWr,
                        i_coeffAdr, o_d, o_dv, nTaps, decFmax)

        inst.convert(hdl=hdl, name=name, path=path)