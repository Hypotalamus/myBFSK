library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use std.textio.all;
use ieee.STD_LOGIC_TEXTIO.all;

use std.env.finish;

library design;
use design.Resampler;

entity tb_Resampler is
end tb_Resampler;

architecture tb of tb_Resampler is

	constant CLK_PERIOD	:	time	:=	10 ns;
    
    signal  i_clk       :   std_logic                       :=  '0';
    signal  i_rst       :   std_logic                       :=  '0';
    signal  i_dv        :   std_logic                       :=  '0';
    signal  i_dPh       :   unsigned(31 downto 0)           :=  (others => '0');
    signal  i_coeffLd   :   std_logic                       :=  '0';
    signal  i_coeff     :   signed (15 downto 0)            :=  (others => '0');
    signal  i_coeffWr   :   std_logic                       :=  '0';
    signal  i_coeffAdr  :   unsigned(15 downto 0)           :=  (others => '0');
    signal  o_dv        :   std_logic                       :=  '0';
    signal  o_next      :   std_logic                       :=  '0';
    signal  i_d_I       :   signed (11 downto 0)            :=  (others => '0');
    signal  i_d_Q       :   signed (11 downto 0)            :=  (others => '0');
    signal  o_d_I       :   signed (11 downto 0)            :=  (others => '0');
    signal  o_d_Q       :   signed (11 downto 0)            :=  (others => '0');

    component Resampler is
        port (
            i_clk: in std_logic;
            i_rst: in std_logic;
            i_dv: in std_logic;
            i_dPh: in unsigned(31 downto 0);
            i_coeffLd: in std_logic;
            i_coeff: in signed (15 downto 0);
            i_coeffWr: in std_logic;
            i_coeffAdr: in unsigned(15 downto 0);
            o_dv: out std_logic;
            o_next: out std_logic;
            i_d_I: in signed (11 downto 0);
            i_d_Q: in signed (11 downto 0);
            o_d_I: out signed (11 downto 0);
            o_d_Q: out signed (11 downto 0)
        );
    end component;

begin

	i_clk	<= not i_clk after CLK_PERIOD/2;

    uut: Resampler
        port map (
            i_clk       =>  i_clk     ,
            i_rst       =>  i_rst     ,
            i_dv        =>  i_dv      ,
            i_dPh       =>  i_dPh     ,
            i_coeffLd   =>  i_coeffLd ,
            i_coeff     =>  i_coeff   ,
            i_coeffWr   =>  i_coeffWr ,
            i_coeffAdr  =>  i_coeffAdr,
            o_dv        =>  o_dv      ,
            o_next      =>  o_next    ,
            i_d_I       =>  i_d_I     ,
            i_d_Q       =>  i_d_Q     ,
            o_d_I       =>  o_d_I     ,
            o_d_Q       =>  o_d_Q     
        );
        
	stimReadIn : process
		file	f			:	text open read_mode is "./stim/iStimuli_Resampler.txt";
		variable    row		:	line;
		variable	stdRd	:	std_logic	:=	'0';
        variable    phRd    :   std_logic_vector(31 downto 0)   :=  (others => '0');
		variable	cRd		:	std_logic_vector(15 downto 0)	:=	(others => '0');
		variable	dRd		:	std_logic_vector(11 downto 0)	:=	(others => '0');
	begin
		wait until falling_edge(i_clk);
		while not endfile(f) loop
			readline(f, row);
			read(row, stdRd);
			i_rst 		    <= stdRd;
			read(row, stdRd);
			i_coeffLd	    <=	stdRd;
			read(row, stdRd);
			i_coeffWr	    <=	stdRd;
			hread(row, cRd);
			i_coeffAdr	    <=	unsigned(cRd);
			hread(row, cRd);
			i_coeff		    <=	signed(cRd);
			hread(row, dRd);
			i_d_I			<=	signed(dRd);
			hread(row, dRd);
			i_d_Q			<=	signed(dRd);
			read(row, stdRd);
			i_dv		    <=	stdRd;
 			hread(row, phRd);
			i_dPh		    <=	unsigned(phRd);           
			wait until rising_edge(i_clk) or falling_edge(i_clk);
		end loop;

        report " **** Test passed **** ";
        finish;
        
	end process;
    
    stimCheck: process
        file f          :   text open read_mode is "./stim/oStimuli_Resampler.txt";
        variable    row :   line;
        variable    dI  :   std_logic_vector(11 downto 0)   :=  (others => '0');
        variable    dQ  :   std_logic_vector(11 downto 0)   :=  (others => '0');
        variable    nxt :   std_logic                       :=  '0';
    begin
        while not endfile(f) loop
            wait until rising_edge(i_clk);
            if o_dv = '1' then
                readline(f, row);
                hread(row, dI);
                hread(row, dQ);
                read(row, nxt);
                
                assert (o_d_I = signed(dI)) and (o_d_Q = signed(dQ)) and (o_next = nxt)
                    report " Expected: (" & to_hstring(dI) & "h, " & to_hstring(dQ) & "h)," &
                           " next = 1'b" & to_string(nxt) & "; " &  
                           " Received: (" & to_hstring(o_d_I) & "h, " & to_hstring(o_d_Q) & "h)," &
                           " o_next = 1'b" & to_string(o_next) & ". "  
                    severity failure;
            end if;
        end loop;
        wait;
        
    end process;
    
end tb;
