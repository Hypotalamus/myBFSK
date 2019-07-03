""" Model of quadrature Numeric Controlled Oscillator """

import numpy as np
from myhdl import *


class NCO(object):
    """ Quadrature NCO model
        List of methods:
        - _quantPh - Quantize phase for integer implementation;
        - _quantOut - Quantize output for integer implementation;
        - setdPhi - Set phase increment at each step;
        - update - Get new NCO sample;
        - rtl - rtl-model of NCO;
        - test_rtl - Testbench. Compare output of rtl-model with 'int' implementation;
        - convert - Convert myHDL-code to HDL file;
     """

    def __init__(self, f0, fs, ampl=1., phi0=0., impl=('float', None)):
        """
            - f0 - nominal output frequency
            - fs - sampling frequency
            - ampl - amplitude of output signal
            - phi0 - initial phase of output signal in [0, 1) (at t = 0)
            - impl - format of data for implementation, in [('float', None), ('int', (phW, ROMadrW, outW))]
                - ('float', None) - process data in floating point format
                - ('int', (phW, ROMadrW, outW)) - process data in floating point format emulating integer format, where
                    - phW - bitwidth of NCO phase, integer > 0
                    - ROMadrW - bitwidth of ROM address with sine/cosine samples (1/4-period ROM)
                    - outW - bitwidth of NCO output sample, integer > 0
                    phi0 and ampl are modified (quantized) accordingly
        """
        self._f0 = f0
        self._fs = fs
        assert impl[0] in ['float', 'int'], (f"Implementation type must be from ['float', 'int']." +
                                             f" Current is '{impl[0]}'")
        self.iType = impl[0]
        if self.iType == 'int':
            self.phW, self.ROMadrW, self.outW = impl[1]
        else:
            self.phW, self.ROMadrW, self.outW = None, None, None
        self.ampl = ampl
        self.phi0 = self._quantPh(phi0)
        self.phi = self._quantPh(phi0)
        self.setdPhi()
        self.rtlRefs = np.array([], dtype=complex)
        self.rtlOuts = np.array([], dtype=complex)

    def reset(self):
        """ Reset NCO to initial state """
        self.phi = self.phi0
        self.rtlRefs = np.array([], dtype=complex)
        self.rtlOuts = np.array([], dtype=complex)

    @property
    def ampl(self):
        return self._ampl

    @ampl.setter
    def ampl(self, val):
        if self.iType == 'int':
            dW = np.ceil(np.log2(val))
            self.outSc = self.outW - dW - 1
        self._ampl = val

    @property
    def f0(self):
        return self._f0

    @f0.setter
    def f0(self, val):
        self._f0 = val
        self.setdPhi()

    @property
    def fs(self):
        return self._fs

    @fs.setter
    def fs(self, val):
        self._fs = val
        self.setdPhi()

    @property
    def dPhi(self):
        return self._dPhi

    def _quantPh(self, ph):
        """ Quantize phase for integer implementation
            - ph - unquantized phase;
        """
        if self.iType == 'float':
            return ph % 1.
        elif self.iType == 'int':
            phQ = np.round(ph * 2 ** self.phW) / 2 ** self.phW
            return phQ % 1.

    def _quantOut(self, val):
        """ Quantize output for integer implementation
            - val - unquantized output value, complex number
            - return complex number without fractional part in real and imaginary components if self.iType == 'int'
        """
        if self.iType == 'float':
            return val
        elif self.iType == 'int':
            valQ = np.round(val * (2 ** self.outSc - 1))
            return valQ

    def setdPhi(self):
        """ Set phase increment at each step """
        self._dPhi = self._quantPh(self.f0 / self.fs)

    def update(self, theta):
        """ Get new NCO sample
            - theta - additional phase
            Return complex number
        """
        thetaQ = self._quantPh(theta)
        ph = np.floor((self.phi + thetaQ) * 2 ** (self.ROMadrW + 2)) / 2 ** (
                    self.ROMadrW + 2) if self.iType == 'int' else \
            self.phi + thetaQ
        val = self.ampl * np.exp(1j * 2 * np.pi * ph)
        self.phi = (self.phi + self.dPhi) % 1.
        return self._quantOut(val)

    @block
    def rtl(self, i_clk, i_rst, i_fw, i_ph0, i_ce, o_cos, o_sin, o_dv, tblDepth):
        """ rtl-model of NCO
            ~~~~ input ports ~~~~
            - i_clk - clock;
            - i_rst - synchronous reset;
            - i_fw - frequency word;
            - i_ph0 - initial phase (after reset), can be set to None;
            - i_ce - clock enable;
            ~~~~ output ports ~~~~
            - o_cos - cosine output;
            - o_sin - sine output;
            - o_dv - output valid, can be set None;
            ~~~~ parameters ~~~~
            tblDepth - Depth of ROM with sine/cosine samples
        """

        acc = Signal(i_fw.val)
        ph = Signal(i_fw.val)
        adr = Signal(intbv(0, min=0, max=tblDepth))
        MAX = adr.max - 1
        MAX_ROM = 2 ** (self.outW - 1)
        ROM1 = tuple(int(self._quantOut(x)) for x in self.ampl * np.cos(np.pi / 2 / tblDepth * np.arange(tblDepth)))
        ROM2 = tuple(int(self._quantOut(x)) for x in self.ampl * np.sin(np.pi / 2 / tblDepth * np.arange(tblDepth)))
        cos, sin = [Signal(intbv(0, min=0, max=MAX_ROM)) for _ in range(2)]
        outCos, outSin = [Signal(o_cos.val) for _ in range(2)]
        outSel, outSelReg = [Signal(intbv(0)[2:]) for _ in range(2)]

        @always(i_clk.posedge)
        def ROMdesc():
            cos.next = ROM1[int(adr)]
            sin.next = ROM2[int(adr)]

        @always_seq(i_clk.posedge, reset=i_rst)
        def intRegs():
            if i_ce:
                acc.next = acc + i_fw
                outSel.next = ph[:len(ph) - 2]

        if i_ph0 is not None:
            @always_seq(i_clk.posedge, reset=i_rst)
            def regs_opt():
                if i_ce:
                    ph.next = acc + i_ph0
        elif i_ph0 is None:
            @always_comb
            def wires_opt():
                ph.next = acc

        @always_comb
        def wires():
            adr.next = ph[len(ph) - 2:len(ph) - 2 - len(adr)]

        @always_seq(i_clk.posedge, reset=i_rst)
        def quadSel():
            if i_ce:
                if outSel == 0:
                    outCos.next = cos
                    outSin.next = sin
                elif outSel == 1:
                    outCos.next = 0 - sin
                    outSin.next = cos
                elif outSel == 2:
                    outCos.next = 0 - cos
                    outSin.next = 0 - sin
                elif outSel == 3:
                    outCos.next = sin
                    outSin.next = 0 - cos

        @always_comb
        def output():
            o_cos.next = outCos
            o_sin.next = outSin

        if o_dv is not None:
            @always_comb
            def out_opt():
                o_dv.next = i_ce

        return instances()

    @block
    def test_rtl(self, period=10, wrToFile=False):
        """ Testbench. Compare output of rtl-model with 'int' implementation
            - period - period of clock, integer;
            - wrToFile - if True, then stimuli and output write to file
        """
        assert self.iType == 'int', (" Testbench of rtl-model is valid only for 'int' implementation. " +
                                     " But current implementation is %s " % self.iType)
        delayRTL = 3
        i_clk = Signal(bool(0))
        i_rst = ResetSignal(0, active=bool(1), isasync=False)
        i_fw = Signal(modbv(0)[self.phW:])
        i_ph0 = Signal(modbv(int(self.phi0 * 2 ** self.phW))[self.phW:])
        i_ce = Signal(bool(0))
        o_cos = Signal(intbv(0, min=-2 ** (self.outW - 1), max=2 ** (self.outW - 1)))
        o_sin = Signal(intbv(0, min=-2 ** (self.outW - 1), max=2 ** (self.outW - 1)))
        o_dv = Signal(bool(0))
        tblDepth = 2 ** self.ROMadrW

        uut = self.rtl(i_clk, i_rst, i_fw, i_ph0, i_ce, o_cos, o_sin, o_dv, tblDepth)

        @always(delay(int(period // 2)))
        def clk_gen():
            i_clk.next = not i_clk

        if wrToFile:
            inFile = open('./vhdl/stim/iStimuli_NCO.txt', 'w')

            @instance
            def wrInStimuli():
                while True:
                    yield i_clk.posedge
                    print(f"{int(i_rst.val)} {i_fw.val} {i_ph0.val} {int(i_ce.val)} ", file=inFile)
                    yield i_clk.negedge
                    print(f"{int(i_rst.val)} {i_fw.val} {i_ph0.val} {int(i_ce.val)} ", file=inFile)

            outFile = open('./vhdl/stim/oStimuli_NCO.txt', 'w')

            @always(i_clk.posedge)
            def wrOut():
                if o_dv:
                    print(f"{o_cos.val} {o_sin.val}", file=outFile)

        @instance
        def stimuli():
            yield i_clk.posedge
            i_rst.next = bool(1)
            for _ in range(3):
                yield i_clk.posedge
            yield i_clk.negedge
            i_rst.next = bool(0)
            i_fw.next = int(self.dPhi * 2 ** self.phW)
            i_ph0.next = int(self.phi0 * 2 ** self.phW)
            i_ce.next = bool(1)
            nTick = 0
            refPnt = 0
            self.rtlRefs = np.array([], dtype=complex)
            self.rtlOuts = np.array([], dtype=complex)
            self.reset()
            for ii in range(2000):
                ref = self.update(0.)
                self.rtlRefs = np.append(self.rtlRefs, ref)
                yield i_clk.negedge
                if o_dv:
                    rtlOut = int(o_cos.val) + 1j * int(o_sin.val)
                    self.rtlOuts = np.append(self.rtlOuts, rtlOut)
                    if nTick < delayRTL - 1:
                        nTick += 1
                    else:
                        assert rtlOut == self.rtlRefs[refPnt], (" Outputs from rtl-model and reference " +
                                                                f"doesn't match: {rtlOut} != {self.rtlRefs[refPnt]}" +
                                                                f", ref pointer = {refPnt}")
                        refPnt += 1
            print(' **** Test passed **** ')
            raise StopSimulation

        return instances()

    def convert(self, isPhi0, isOdv, hdl='VHDL', name='NCO', path='./vhdl'):
        """ Convert myHDL-code to HDL file
            - isPhi0 - add (True) or not (False) i_phi0 input, boolean
            - isOdv - add (True) or not (False) o_dv output, boolean
            - hdl - HDL language in ['VHDL', 'Verilog']
            - name - name of the entity (module) and appropriate file, string
            - path - destination folder, string
        """
        i_clk = Signal(bool(0))
        i_rst = ResetSignal(0, active=bool(1), isasync=False)
        i_fw = Signal(intbv(0)[self.phW:])
        if isPhi0:
            i_ph0 = Signal(intbv(0)[self.phW:])
        else:
            i_ph0 = None
        i_ce = Signal(bool(0))
        o_cos = Signal(intbv(0, min=-2 ** (self.outW - 1), max=2 ** (self.outW - 1)))
        o_sin = Signal(intbv(0, min=-2 ** (self.outW - 1), max=2 ** (self.outW - 1)))
        if isOdv:
            o_dv = Signal(bool(0))
        else:
            o_dv = None
        tblDepth = 2 ** self.ROMadrW

        inst = self.rtl(i_clk, i_rst, i_fw, i_ph0, i_ce, o_cos, o_sin, o_dv, tblDepth)
        inst.convert(hdl=hdl, name=name, path=path)