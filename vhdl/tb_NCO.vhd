library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use std.textio.all;
use ieee.std_logic_textio.all;

use std.env.finish;

library design;
use design.NCO;

entity tb_NCO is
end tb_NCO;

architecture tb of tb_NCO is

    constant CLK_PERIOD :   time    :=  10 ns;
    
    signal  i_clk   :   std_logic                       :=  '0';
    signal  i_rst   :   std_logic                       :=  '0';
    signal  i_fw    :   std_logic_vector(31 downto 0)   :=  (others => '0');
    signal  i_ph0   :   std_logic_vector(31 downto 0)   :=  (others => '0');
    signal  i_ce    :   std_logic                       :=  '0';
    signal  o_cos   :   signed(11 downto 0)             :=  (others => '0');
    signal  o_sin   :   signed(11 downto 0)             :=  (others => '0');
    signal  o_dv    :   std_logic                       :=  '0';

    component NCO is
        port (
            i_clk: in std_logic;
            i_rst: in std_logic;
            i_fw: in unsigned(31 downto 0);
            i_ph0: in unsigned(31 downto 0);
            i_ce: in std_logic;
            o_cos: out signed (11 downto 0);
            o_sin: out signed (11 downto 0);
            o_dv: out std_logic
        );
    end component;
    
begin

    i_clk   <=  not i_clk   after CLK_PERIOD/2;

    uut: NCO
        port map (
            i_clk   =>  i_clk,
            i_rst   =>  i_rst,
            i_fw    =>  unsigned(i_fw),
            i_ph0   =>  unsigned(i_ph0),
            i_ce    =>  i_ce ,
            o_cos   =>  o_cos,
            o_sin   =>  o_sin,
            o_dv    =>  o_dv
        );

    stimReadIn: process
        file    f               :   text open read_mode is "./stim/iStimuli_NCO.txt";
        variable    row         :   line;
        variable    stdRd       :   std_logic                       :=  '0';
        variable    stdv32Rd    :   std_logic_vector(31 downto 0)   :=  (others => '0');
    begin
        wait until falling_edge(i_clk);
        while not endfile(f) loop
            readline(f, row);
            read(row, stdRd);
            i_rst   <=  stdRd;
            hread(row, stdv32Rd);
            i_fw    <=  stdv32Rd;
            hread(row, stdv32Rd);
            i_ph0   <=  stdv32Rd;
            read(row, stdRd);
            i_ce    <=  stdRd;
            wait until rising_edge(i_clk) or falling_edge(i_clk);
        end loop;
        
        report " **** Test passed **** ";
        finish;
        
    end process;
    
    stimCheck: process
        file    f       :   text open read_mode is "./stim/oStimuli_NCO.txt";
        variable    row :   line;
        variable    cos :   std_logic_vector(11 downto 0)   :=  (others =>  '0');
        variable    sin :   std_logic_vector(11 downto 0)   :=  (others =>  '0');
    begin
        while not endfile(f) loop
            wait until rising_edge(i_clk);
            if o_dv = '1' then
                readline(f, row);
                hread(row, cos);
                hread(row, sin);
                
                assert (o_cos = signed(cos)) and (o_sin = signed(sin))
                    report "Expected: (" & to_hstring(cos) & "h, " & to_hstring(sin) & "h), " &
                           "Received: (" & to_hstring(o_cos) & "h, " & to_hstring(o_sin) & "h)."
                    severity failure;                
            end if;
        end loop;
        wait;
        
    end process;
        
end tb;