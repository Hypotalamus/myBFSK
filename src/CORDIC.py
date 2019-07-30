""" Model of CORDIC for angle calculation """

import numpy as np
from myhdl import *
from .utils import IQdata

class CORDIC(object):
    """ Processor for complex number argument calculation via CORDIC-algorithm
        List of methods:
        - setAngs - Calculate angles for each stage;
        - calc - Calculate angle of pnt;
        - norm - Normalize number to improve precision in angle calculation
            (myHDL);
        - rotate - One stage of rotate operation (myHDL);
        - rtl - RTL-model of CORDIC in myHDL;
        - test_rtl - Testbench. Compare output of rtl-model with 'int' implementation;
        - convert - Convert myHDL-code to HDL file;
     """

    def __init__(self, dBw, cordBw, stQnt, cordBwInt):
        """
            - dBw - input data bitwidth
            - cordBw - bitwidth of output angle
            - stQnt - number of CORDIC stages
            - cordBwInt - internal angle bitwidth >= self.cordBw, if None then = self.cordBw
        """
        self.dBw = dBw
        self.cordBw = cordBw
        self.stQnt = stQnt
        if cordBwInt:
            self.cordBwInt = cordBwInt
        else:
            self.cordBwInt = self.cordBw
        self.setAngs()
        self.anglesRef = np.array([], dtype=int)
        self.rtlOuts = np.array([], dtype=int)

    def setAngs(self, stQnt=None, cordBwInt=None):
        """ Calculate angles for each stage
            - stQnt - number of CORDIC stages, if None, old value remain
            - cordBwInt - internal angle bitwidth, if None, old value remain
        """
        if stQnt:
            self.stQnt = stQnt
        if cordBwInt:
            self.cordBwInt = cordBwInt
        self.angs = np.array([int(np.round(2 ** (self.cordBwInt - 1) / np.pi * np.arctan(1 / 2 ** ii)))
                              for ii in range(self.stQnt)])

    def calc(self, pnt):
        """ Calculate angle of pnt
            - pnt - complex number with integer (without fractional part) real and imaginary components
            return angle - integer number in [-2**(cordBw-1), 2**(cordBw-1)) where 2**(cordBw-1) corresponds to np.pi
        """
        pntInt = pnt
        x, y = int(pntInt.real), int(pntInt.imag)
        while (-2 ** (self.dBw - 2) <= x < 2 ** (self.dBw - 2)) and (-2 ** (self.dBw - 2) <= y < 2 ** (self.dBw - 2)):
            x <<= 1
            y <<= 1
        x <<= self.cordBwInt - self.dBw
        y <<= self.cordBwInt - self.dBw
        xprev = x
        z = 0
        # initial rotate by +/- np.pi/2
        if x >= 0:
            z = 0
        elif y >= 0:
            x = y
            y = -xprev
            xprev = x
            z = 2 ** (self.cordBwInt - 2)
        else:
            x = -y
            y = xprev
            xprev = x
            z = - 2 ** (self.cordBwInt - 2)

        for ii in range(self.stQnt):
            if y >= 0:
                x = x + (y >> ii)
                y = y - (xprev >> ii)
                xprev = x
                z += self.angs[ii]
            elif y < 0:
                x = x - (y >> ii)
                y = y + (xprev >> ii)
                xprev = x
                z -= self.angs[ii]
        if self.cordBwInt > self.cordBw:
            tmp = (z + (1 << self.cordBwInt - self.cordBw - 1)) >> (self.cordBwInt - self.cordBw)
            if tmp >= 2 ** (self.cordBw - 1):
                return 2 ** (self.cordBw - 1) - 1
            else:
                return tmp
        else:
            return z

    @block
    def norm(self, i_clk, i_rst, i_d, i_dv, o_d, o_dv):
        """ Normalize number to improve precision in angle calculation
            ~~~~ input ports ~~~~
            - i_clk - clock;
            - i_rst - synchronous reset;
            - i_d - input complex data, class IQdata();
            - i_dv - input data valid flag;
            ~~~~ output ports ~~~~
            - o_d - output data, class IQData();
            - o_dv - output data valid
        """
        dRe, dIm, qRe, qIm = [Signal(i_d.I.val) for _ in range(4)]
        reSig, imSig, reSSB, imSSB, val = [Signal(bool(0)) for _ in range(5)]

        @always_comb
        def getInput():
            dRe.next = i_d.I
            dIm.next = i_d.Q
            reSig.next = i_d.I[len(i_d.I) - 1]
            imSig.next = i_d.Q[len(i_d.Q) - 1]
            reSSB.next = i_d.I[len(i_d.I) - 2]
            imSSB.next = i_d.Q[len(i_d.Q) - 2]

        @always_seq(i_clk.posedge, i_rst)
        def normLogic():
            if i_dv:
                val.next = bool(1)
                if (reSig == reSSB) and (imSig == imSSB):
                    qRe.next = dRe << 1
                    qIm.next = dIm << 1
                else:
                    qRe.next = dRe
                    qIm.next = dIm

        @always_comb
        def output():
            o_d.I.next = qRe
            o_d.Q.next = qIm
            o_dv.next = val and i_dv

        return instances()

    @block
    def rotate(self, i_clk, i_rst, i_d, i_dv, i_ang, o_d, o_dv, o_ang, dAng, stInd):
        """ One stage of rotate operation
            ~~~~ input ports ~~~~
            - i_clk - clock;
            - i_rst - synchronous reset;
            - i_d - input data, class IQdata();
            - i_dv - input data valid flag;
            - i_ang - input angle;
            ~~~~ output ports ~~~~
            - o_d - output data, class IQdata();
            - o_dv - output data valid flag;
            - o_ang - output angle;
            ~~~~ parameters ~~~~
            - dAng - rotational angle, natural
            - stInd - stage index, natural
        """
        dRe, dIm, qRe, qIm = [Signal(i_d.I.val) for _ in range(4)]
        ang, qAng = [Signal(i_ang.val) for _ in range(2)]
        val = Signal(bool(0))

        @always_comb
        def getInput():
            dRe.next = i_d.I
            dIm.next = i_d.Q
            ang.next = i_ang

        @always_seq(i_clk.posedge, i_rst)
        def rotLogic():
            if i_dv:
                val.next = bool(1)
                if dIm >= 0:
                    qRe.next = dRe + (dIm >> stInd)
                    qIm.next = dIm - (dRe >> stInd)
                    qAng.next = ang + dAng
                else:
                    qRe.next = dRe - (dIm >> stInd)
                    qIm.next = dIm + (dRe >> stInd)
                    qAng.next = ang - dAng

        @always_comb
        def output():
            o_d.I.next = qRe
            o_d.Q.next = qIm
            o_ang.next = qAng
            o_dv.next = val and i_dv

        return instances()

    @block
    def rtl(self, i_clk, i_rst, i_d, i_dv, o_ang, o_dv):
        """ RTL-model of CORDIC in myHDL
            ~~~~ input ports ~~~~
            - i_clk - clock;
            - i_rst - synchronous reset;
            - i_d - input data, class IQdata();
            - i_dv - input data valid flag;
            ~~~~ output ports ~~~~
            - o_ang - output angle;
            - o_dv - output angle valid flag
        """
        dRe, dIm = [Signal(i_d.I.val) for _ in range(2)]

        @always_comb
        def getInputs():
            dRe.next = i_d.I
            dIm.next = i_d.Q

        # normalize pipeline
        normBlocks = []
        for ii in range(self.dBw - 2):
            normOut = IQdata(self.dBw)
            normOutV = Signal(bool(0))
            if ii == 0:
                normIn = IQdata(self.dBw)
                normIn.I = dRe
                normIn.Q = dIm
                normBlocks.append(self.norm(i_clk, i_rst, normIn, i_dv, normOut, normOutV))
                normIn = normOut
                normInV = normOutV
            else:
                normBlocks.append(self.norm(i_clk, i_rst, normIn, normInV, normOut, normOutV))
                normIn = normOut
                normInV = normOutV

        # expand normalized data to internal bitwidth
        dReExp, dImExp = [Signal(intbv(0, min=-2 ** (self.cordBwInt + 1), max=2 ** (self.cordBwInt + 1))) for _ in
                          range(2)]
        dExpVal = Signal(bool(0))
        angInit = Signal(intbv(0, min=-2 ** (self.cordBwInt), max=2 ** (self.cordBwInt)))

        @always_comb
        def expand():
            dReExp.next = normIn.I << (self.cordBwInt - self.dBw)
            dImExp.next = normIn.Q << (self.cordBwInt - self.dBw)
            dExpVal.next = normInV

        # initial rotate by +/- np.pi/2
        fRotD = IQdata(self.cordBwInt + 2)
        fRotVdel = Signal(bool(0))
        fRotV = Signal(bool(0))

        @always_seq(i_clk.posedge, i_rst)
        def firstRot():
            if dExpVal:
                fRotVdel.next = 1
                if dReExp >= 0:
                    fRotD.I.next = dReExp
                    fRotD.Q.next = dImExp
                    angInit.next = 0
                elif dImExp >= 0:
                    fRotD.I.next = dImExp
                    fRotD.Q.next = -dReExp
                    angInit.next = 2 ** (self.cordBwInt - 2)
                else:
                    fRotD.I.next = -dImExp
                    fRotD.Q.next = dReExp
                    angInit.next = - 2 ** (self.cordBwInt - 2)

        @always_comb
        def setfRotV():
            fRotV.next = fRotVdel and dExpVal

        # Another rotates
        rotBlocks = []
        for ii, rAng in enumerate(self.angs):
            rotOut = IQdata(self.cordBwInt + 2)
            rotAngO = Signal(intbv(0, min=-2 ** (self.cordBwInt), max=2 ** (self.cordBwInt)))
            rotOutV = Signal(bool(0))
            if ii == 0:
                rotBlocks.append(self.rotate(i_clk, i_rst, fRotD, fRotV, angInit, rotOut, rotOutV, rotAngO,
                                             int(rAng), ii))
                rotIn = rotOut
                rotAngI = rotAngO
                rotInV = rotOutV
            else:
                rotBlocks.append(self.rotate(i_clk, i_rst, rotIn, rotInV, rotAngI, rotOut, rotOutV, rotAngO,
                                             int(rAng), ii))
                rotIn = rotOut
                rotAngI = rotAngO
                rotInV = rotOutV

                # round output
        if self.cordBwInt > self.cordBw:
            p5 = 1 << self.cordBwInt - self.cordBw - 1
            angNorm = Signal(intbv(0, min=rotAngI.min, max=rotAngI.max))

            @always_comb
            def getAngNorm():
                angNorm.next = (rotAngI + p5) >> (self.cordBwInt - self.cordBw)

            @always_comb
            def outWRound():
                if angNorm == 2 ** (self.cordBw - 1):
                    o_ang.next = 2 ** (self.cordBw - 1) - 1
                else:
                    o_ang.next = angNorm
                o_dv.next = rotInV

        else:
            @always_comb
            def outWoRound():
                o_ang.next = rotAngI
                o_dv.next = rotInV

        return instances()

    @block
    def test_rtl(self, period=10, wrToFile=False, randQnt=1000):
        """ Testbench. Compare output of rtl-model with 'int' implementation
            - period - period of clock, integer;
            - wrToFile - if True, then stimuli and output write to file
            - randQnt - number of test input data
        """
        i_clk = Signal(bool(0))
        i_rst = ResetSignal(0, active=bool(1), isasync=False)
        i_d = IQdata(self.dBw)
        i_dv = Signal(bool(0))
        o_ang = Signal(intbv(0, min=-2 ** (self.cordBw - 1), max=2 ** (self.cordBw - 1)))
        o_dv = Signal(bool(0))

        uut = self.rtl(i_clk, i_rst, i_d, i_dv, o_ang, o_dv)

        @always(delay(int(period // 2)))
        def clk_gen():
            i_clk.next = not i_clk

        if wrToFile:
            inFile = open('./vhdl/stim/iStimuli_CORDIC.txt', 'w')

            @instance
            def wrInStimuli():
                while True:
                    yield i_clk.posedge
                    print(f"{int(i_rst.val)} {i_d.I.val} {i_d.Q.val} {int(i_dv.val)} ", file=inFile)
                    yield i_clk.negedge
                    print(f"{int(i_rst.val)} {i_d.I.val} {i_d.Q.val} {int(i_dv.val)} ", file=inFile)

            outFile = open('./vhdl/stim/oStimuli_CORDIC.txt', 'w')

            @always(i_clk.posedge)
            def wrOut():
                if o_dv:
                    print(f"{o_ang.val} {int(o_dv.val)}", file=outFile)

        compMax = 2 ** (self.dBw - 1) - 1
        pnts = np.random.choice([-1, 1], randQnt) * np.random.randint(0, compMax, size=randQnt) + \
               1j * np.random.choice([-1, 1], randQnt) * np.random.randint(0, compMax, size=randQnt)
        self.anglesRef = np.array([self.calc(pnt) for pnt in pnts]).astype(int)
        self.rtlOuts = np.array([], dtype=int)

        @instance
        def stimuli():
            yield i_clk.posedge
            i_rst.next = bool(1)
            for _ in range(3):
                yield i_clk.posedge
            yield i_clk.negedge
            i_rst.next = bool(0)

            i_dv.next = bool(1)
            for pnt in pnts:
                i_d.I.next = int(pnt.real)
                i_d.Q.next = int(pnt.imag)
                yield i_clk.negedge
                if o_dv:
                    rtlOut = int(o_ang.val)
                    self.rtlOuts = np.append(self.rtlOuts, rtlOut)

            assert np.all(self.rtlOuts == self.anglesRef[:len(self.rtlOuts)]), (
                " Outputs from rtl-model and reference doesn't match ")

            print(' **** Test passed **** ')
            raise StopSimulation

        return instances()

    def convert(self, hdl='VHDL', name='CORDIC', path='./vhdl'):
        """ Convert myHDL-code to HDL file
            - hdl - HDL language in ['VHDL', 'Verilog']
            - name - name of the entity (module) and appropriate file, string
            - path - destination folder, string
        """
        i_clk = Signal(bool(0))
        i_rst = ResetSignal(0, active=bool(1), isasync=False)
        i_d = IQdata(self.dBw)
        i_dv = Signal(bool(0))
        o_ang = Signal(intbv(0, min=-2 ** (self.cordBw - 1), max=2 ** (self.cordBw - 1)))
        o_dv = Signal(bool(0))

        inst = self.rtl(i_clk, i_rst, i_d, i_dv, o_ang, o_dv)
        inst.convert(hdl=hdl, name=name, path=path)