library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use std.textio.all;
use ieee.STD_LOGIC_TEXTIO.all;

use std.env.finish;

entity tb_PolyphaseDecimator is
end tb_PolyphaseDecimator;

architecture tb of tb_PolyphaseDecimator is

	constant CLK_PERIOD	:	time	:=	10 ns;

	signal	i_clk		:	std_logic						:=	'0';
	signal	i_rst		:	std_logic						:=	'0';
	signal	i_coeffLd	:	std_logic						:=	'0';
	signal	i_coeffWr	:	std_logic						:=	'0';
	signal	i_coeffAdr	:	std_logic_vector(7 downto 0)	:=	(others => '0');
	signal	i_coeff		:	std_logic_vector(15 downto 0)	:=	(others => '0');
	signal	i_d_I		:	std_logic_vector(11 downto 0)	:=	(others => '0');
	signal	i_d_Q		:	std_logic_vector(11 downto 0)	:=	(others => '0');
	signal	i_dv		:	std_logic						:=	'0';
    signal  i_decF      :   std_logic_vector(3 downto 0)    :=  (others => '0');
    signal  o_d_I       :   signed(11 downto 0)             :=  (others => '0');
    signal  o_d_Q       :   signed(11 downto 0)             :=  (others => '0');
    signal  o_dv        :   std_logic                       :=  '0';
    
    component PolyphaseDecimator is
        port (
            i_clk: in std_logic;
            i_rst: in std_logic;
            i_dv: in std_logic;
            i_decF: in unsigned(3 downto 0);
            i_coeffLd: in std_logic;
            i_coeff: in signed (15 downto 0);
            i_coeffWr: in std_logic;
            i_coeffAdr: in unsigned(7 downto 0);
            o_dv: out std_logic;
            i_d_I: in signed (11 downto 0);
            i_d_Q: in signed (11 downto 0);
            o_d_I: out signed (11 downto 0);
            o_d_Q: out signed (11 downto 0)
        );
    end component;

begin

	i_clk	<= not i_clk after CLK_PERIOD/2;

    uut: PolyphaseDecimator
        port map (
            i_clk       =>  i_clk,
            i_rst       =>  i_rst,
            i_d_I       =>  signed(i_d_I),
            i_d_Q       =>  signed(i_d_Q),
            i_dv        =>  i_dv,
            i_decF      =>  unsigned(i_decF),
            i_coeffLd   =>  i_coeffLd,
            i_coeff     =>  signed(i_coeff),
            i_coeffWr   =>  i_coeffWr,
            i_coeffAdr  =>  unsigned(i_coeffAdr),
            o_d_I       =>  o_d_I,
            o_d_Q       =>  o_d_Q,    
            o_dv        =>  o_dv
        );


	stimReadIn : process
		file	f			:	text open read_mode is "./stim/iStimuli_Decimator.txt";
		variable    row		:	line;
		variable	stdRd	:	std_logic	:=	'0';
		variable	cAdrRd	:	std_logic_vector(7 downto 0)	:=	(others => '0');
		variable	cRd		:	std_logic_vector(15 downto 0)	:=	(others => '0');
		variable	dRd		:	std_logic_vector(11 downto 0)	:=	(others => '0');
        variable    dfRd    :   std_logic_vector(3 downto 0)    :=  (others => '0');
	begin
		wait until falling_edge(i_clk);
		while not endfile(f) loop
			readline(f, row);
			read(row, stdRd);
			i_rst 		<= stdRd;
			read(row, stdRd);
			i_coeffLd	<=	stdRd;
			read(row, stdRd);
			i_coeffWr	<=	stdRd;
			hread(row, cAdrRd);
			i_coeffAdr	<=	cAdrRd;
			hread(row, cRd);
			i_coeff		<=	cRd;
			hread(row, dRd);
			i_d_I			<=	dRd;
			hread(row, dRd);
			i_d_Q			<=	dRd;
			read(row, stdRd);
			i_dv		<=	stdRd;
 			hread(row, dfRd);
			i_decF		<=	dfRd;           
			wait until rising_edge(i_clk) or falling_edge(i_clk);
		end loop;

        report " **** Test passed **** ";
        finish;
        
	end process;
    
    stimCheck: process
        file f          :   text open read_mode is "./stim/oStimuli_Decimator.txt";
        variable    row :   line;
        variable    dI  :   std_logic_vector(11 downto 0)   :=  (others => '0');
        variable    dQ  :   std_logic_vector(11 downto 0)   :=  (others => '0');
    begin
        while not endfile(f) loop
            wait until rising_edge(i_clk);
            if o_dv = '1' then
                readline(f, row);
                hread(row, dI);
                hread(row, dQ);
                
                assert (o_d_I = signed(dI)) and (o_d_Q = signed(dQ))
                    report " Expected: (" & to_hstring(dI) & "h, " & to_hstring(dQ) & "h)," &
                           " Received: (" & to_hstring(o_d_I) & "h, " & to_hstring(o_d_Q) & "h)."
                    severity failure;
            end if;
        end loop;
        wait;
        
    end process;
    
end tb;
