library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use std.textio.all;
use ieee.std_logic_textio.all;

use std.env.finish;

library design;
use design.CORDIC;

entity tb_CORDIC is
end tb_CORDIC;

architecture tb of tb_CORDIC is

    constant CLK_PERIOD :   time    :=  10 ns;
    
    signal  i_clk   :   std_logic                       :=  '0';
    signal  i_rst   :   std_logic                       :=  '0';
    signal  i_d_I   :   std_logic_vector(11 downto 0)   :=  (others => '0');
    signal  i_d_Q   :   std_logic_vector(11 downto 0)   :=  (others => '0');
    signal  i_dv    :   std_logic                       :=  '0';
    signal  o_ang   :   signed(15 downto 0)             :=  (others => '0');
    signal  o_dv    :   std_logic                       :=  '0';

    component CORDIC is
        port (
            i_clk: in std_logic;
            i_rst: in std_logic;
            i_dv: in std_logic;
            o_ang: out signed (15 downto 0);
            o_dv: out std_logic;
            i_d_I: in signed (11 downto 0);
            i_d_Q: in signed (11 downto 0)
        );
    end component;
    
begin

    i_clk   <=  not i_clk   after CLK_PERIOD/2;

    uut: CORDIC
        port map (
            i_clk   =>  i_clk,
            i_rst   =>  i_rst,
            i_dv    =>  i_dv,
            o_ang   =>  o_ang,
            o_dv    =>  o_dv,
            i_d_I   =>  signed(i_d_I),
            i_d_Q   =>  signed(i_d_Q)
        );

    stimReadIn: process
        file    f               :   text open read_mode is "./stim/iStimuli_CORDIC.txt";
        variable    row         :   line;
        variable    stdRd       :   std_logic                       :=  '0';
        variable    stdv12Rd    :   std_logic_vector(11 downto 0)   :=  (others => '0');
    begin
        wait until falling_edge(i_clk);
        while not endfile(f) loop
            readline(f, row);
            read(row, stdRd);
            i_rst   <=  stdRd;
            hread(row, stdv12Rd);
            i_d_I   <=  stdv12Rd;
            hread(row, stdv12Rd);
            i_d_Q   <=  stdv12Rd;
            read(row, stdRd);
            i_dv    <=  stdRd;
            wait until rising_edge(i_clk) or falling_edge(i_clk);
        end loop;
        
        report " **** Test passed **** ";
        finish;
        
    end process;
    
    stimCheck: process
        file    f       :   text open read_mode is "./stim/oStimuli_CORDIC.txt";
        variable    row :   line;
        variable    ang :   std_logic_vector(15 downto 0)   :=  (others =>  '0');
    begin
        while not endfile(f) loop
            wait until rising_edge(i_clk);
            if o_dv = '1' then
                readline(f, row);
                hread(row, ang);
                
                assert (o_ang = signed(ang))
                    report "Expected: (" & to_hstring(ang) & "h, " &
                           "Received: (" & to_hstring(o_ang) & "h)."
                    severity failure;                
            end if;
        end loop;
        wait;
        
    end process;
        
end tb;