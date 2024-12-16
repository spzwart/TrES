# TRES
TRiple Evolution Simulator 

### Description
TRES is a numerical framework for simulating hierarchical triple systems with stellar and planetary components. 
Mass transfer from one star to another and the consequential effect to the orbital dynamics is realized via heuristic recipes.
These recipes are combined with  three-body  dynamics and stellar evolution inluding their mutual influences. 

TRES includes the effects of common-envelope evolution, circularized stable mass transfer, tides, gravitational wave emission and up-to-date stellar evolution through SeBa. Other stellar evolution codes such as SSE or MESA can also be used. Coming soon: transition to N-body calculations (including stellar evolution and dissipative processes) when the system's evolution is not secular anymore. 

This document contains the following parts:

[Compilation](#Compilation)

[Simple examples](#Simple-examples-of-runs)

[TRES-Exo for exoplanet research](#TRES-Exo-for-exoplanet-research)

[TRES with MESA](#TRES-with-MESA)

[Understanding the TRES output](#Understanding-the-TRES-output)

[Reducing the TRES output](#Reducing-the-TRES-output)

[TRES development team](#TRES-development-team)

[References](#References)



## Compilation

TRES makes use of the Astrophysical Multipurpose Software Environment (AMUSE) See https://amusecode.github.io/ for how to install AMUSE. 
Note that for standard TRES usage, the only necessary community code to install is SeBa. 

Thus, after installing the AMUSE pre-requisites, we can simply install the minimal framework and then add SeBa:

```

pip install [--user] amuse-framework
pip install [--user] amuse-<seba>

```

After compiling AMUSE, TRES needs to be installed and compiled by means of the Makefile as following:

First, clone the TRES github repository:

```

git clone https://github.com/amusecode/TRES.git

```

Then, from the root of the cloned respository compile the Makefile:

```

cd seculartriple_TPS
make clean
make 

```

Note: If the above line doesn't work, which may be the case for older versions of amuse, try:
```
(cd seculartriple_TPS)
mv Makefile Makefile_new
mv Makefile_old Makefile
make clean
make 

```



## Simple examples of runs

To evolve a single system with the parameters:
primary mass M=1.2 Solar mass, 
secondary mass m=0.5 Solar mass, 
tertiary mass m=0.6 Solar mass, 
inner and outer eccentricity e=0.1 & 0.5, 
inner and outer orbital separation a=200 & 20000 Solar radii, 
metallicity z=0.001, 
and time T=10 Myrs,  you need to run:

```

python TRES.py -M 1.2 -m 0.5 -l 0.6 -E 0.1 -e 0.5 -A 200 -a 20000 -z 0.001 -T 10 

```
assuming AMUSE is loaded in your python. 

### Input parameters 

The full list of possible input parameters is:
Depreciated (yet still functioning) parameters are given in {}. 

```
                  parameter                               unit / default
--M1    {-M}      inner_primary_mass                      in Solar mass
--M2    {-m}      inner secondary mass                    in Solar mass 
--M3    {-l}      outer mass                              in Solar mass
--Ain   {-A}      inner semi major axis                   in Solar radius
--Aout  {-a}      outer semi major axis                   in Solar radius
--Ein   {-E}      inner eccentricity 
--Eout  {-e}      outer eccentricity 
-i, -I            relative inclination                    in rad  
--Gin   {-G}      inner argument of pericenter            in rad
--Gout  {-g}      outer argument of pericenter            in rad
--Oin   {-O}      inner longitude of ascending node       in rad
                  (outer longitude of ascending nodes = inner - pi)               
-Z      {-z}      metallicity                             default = 0.02 (Solar)
-t, -T            end time                                in Myr
-N, -n            integer number asigned to the triple    default = 0
    
    
-f                name of output file                     default = TRES
-F                type of output file (hdf5/txt)          default = hdf5
--dir_plots       directory for plots for debugging default = "" (current directory)
                  mode  (aka REPORT_DEBUG == True)

--CE    which type of modelling for common envelope evolution   default = 2
        options:
        0:  alpha-ce + alpha-dce
        1:  gamma-ce + alpha-dce
        2:  seba style; combination of gamma-ce, alpha-ce & alpha-dce  
Note that in all cases a double common-envelope is calculated using the alpha-ce.           

--SN_kick_distr   supernova kick distribution   default = 10
        options:
        0:  No kick 
        1:  Hobbs, Lorimer, Lyne & Kramer 2005, 360, 974  
        2:  Hobbs, Lorimer, Lyne & Kramer 2005, 360, 974  scaled down for bh by mass
        3:  Arzoumanian ea 2002, 568, 289
        4:  Arzoumanian ea 2002, 568, 289 scaled down for bh by mass
        5:  Hansen & Phinney 1997, 291, 569
        6:  Hansen & Phinney 1997, 291, 569 scaled down for bh by mass
        7:  Paczynski 1990, 348, 485
        8:  Paczynski 1990, 348, 485 scaled down for bh by mass
        9:  Verbunt, Igoshev & Cator, 2017, 608, 57
        10:  Verbunt, Igoshev & Cator, 2017, 608, 57 scaled down for bh by mass 
        
--max_CPU_time    maximum CPU time allowed (only works in combination with "stop_at_CPU_time")    
                                                default = 3600 (seconds)

```

Additionally, there is a list of stopping conditions that determines whether the simulation of a system should stop at a certain evolutionary phase. 
By default, these stopping conditions are set to True, which means they are in effect. However, the four specific mass transfer cases (stable, unstable, eccentric stable & eccentric unstable) are set to False by default. Once "--no_stop_at_mass_transfer" is set to False, it is possible to set the specific mass transfer cases to True.

```

action items                                    add these to:
--no_stop_at_mass_transfer                      avoid stopping the simulation at the onset of mass transfer 
--no_stop_at_init_mass_transfer                 avoid stopping the simulation if there is mass transfer initially
--no_stop_at_outer_mass_transfer                avoid stopping the simulation when tertiary initiates mass transfer 
                                                methodology is as of yet non-existent
--stop_at_stable_mass_transfer                  avoid stopping the simulation at the onset of stable mass transfer in the inner binary
--stop_at_unstable_mass_transfer                avoid stopping the simulation at the onset of unstable mass transfer in the inner binary (leading to common-envelope evolution)
--stop_at_eccentric_stable_mass_transfer        avoid stopping the simulation at the onset of stable mass transfer in the inner binary if the orbit is eccentric
--stop_at_eccentric_unstable_mass_transfer      avoid stopping the simulation at the onset of unstable mass transfer in the inner binary if the orbit is eccentric
--no_stop_at_merger                             avoid stopping the simulation after a merger
--no_stop_at_inner_collision                    avoid stopping the simulation after a collision in the inner binary
--no_stop_at_outer_collision                    avoid stopping the simulation after a collision involving the outer star
--no_stop_at_disintegrated                      avoid stopping after the system disintegrated into seperate systems
--stop_at_semisecular_regime                    to stop the simulation if the sytem is in the semi secular regime
--stop_at_SN                                    to stop the simulation when a supernova occurs
--stop_at_CPU_time                              to stop the simulation when the computational time exceeds a given value

```


### Multiple systems with specified parameters

If you need to follow the triple evolution for multiple systems with parameters which are already specified you can start TRES multiple times, e.g.
```

python TRES.py -M 1.2 -m 0.5 -l 0.6 -E 0.1 -e 0.5 -A 200 -a 20000 -z 0.001 -T 10 
python TRES.py -M 1.5 -m 1 -l 0.6 -E 0.1 -e 0.5 -A 50 -a 20000 -z 0.001 -T 10 
python TRES.py -M 1.5 -m 1 -l 0.05 -E 0.1 -e 0.5 -A 50 -a 20000 -z 0.02 -T 10 

```
This is probably not handy for more than 5 systems. Although this can be added in e.g. a shell or Python script.


### Random population
A random population can be generated with TPS.py with a Monte Carlo based approach, e.g. 
```

python TPS.py -n 5 
python TPS.py -n 10 --M_max 5 --M_min 4  --M_distr 0 --A_max 2000 --A_min 200 --A_distr 2

```

The full list of options is [default]:
Depreciated (yet still functioning) parameters are given in {}. 
```

--M1_max       {--M_max}    upper limit for the inner primary mass [100 Msun]
--M1_min       {--M_min}    lower limit for the inner primary mass [0.1 Msun]
--M1_distr     {--M_distr}  mass function option: 
        0: "Kroupa", #default
        1: "Scalo",
        2: "Miller & Scalo",
        3: "Salpeter",
        4: "Logarithmically flat",
        5: "Eggleton",
        6: "Kroupa for massive stars M>0.5 powerlaw with exp=-2.3",
--Qin_max      {--Q_max}    upper limit for the inner mass ratio [1.]
--Qin_min      {--Q_min}    lower limit for the inner mass ratio [0.]
--Qin_distr    {--Q_distr}  inner mass ratio option: 
       0: "Uniform distribution", #default
       1: "Kroupa IMF",
       2: "Galicher et al. 2016 powerlaw (M^-1.31)", #draws from mass distribution instead of mass ratio distribution, 
--Qout_max     {--q_max}    upper limit for the outer mass ratio [1.]
--Qout_min     {--q_min}    lower limit for the mass of the outer star [0.]
--Qout_distr   {--q_distr}  outer mass ratio option: 
       0: "Uniform distribution", #default
       1: "Kroupa IMF",
       2: "Galicher et al. 2016 powerlaw (M^-1.31)", #draws from mass distribution instead of mass ratio distribution, 
--Ain_max      {--A_max}    upper limit for the inner semi-major axis [5e6 RSun]
--Ain_min      {--A_min}    lower limit for the inner semi-major axis [0.5 RSun]
        Note that the true minimum separation is also affected by RLOF. 
        By default contact or semi-detached systems on the ZAMS are removed from the sample. 
        We recommend to add the command line option --include_circularisation_during_preMS, such that the eccentricity is redrawn (<10x). 
--Ain_distr    {--A_distr}  inner semi-major axis option: 
        0: "Log Uniform distribution", #default
        1: "Constant semi-major axis",
        2: "Tokovinin lognormal mu = 10^5d, sigma = 2.3",
        3: "Lognormal mu = 10^3.5d, sigma = 2.3",
        4: "Rizzuto Lognormal mu = 10^0.95 AU, sigma = 1.35",
        5: "Sana et al. 2012",
        6: "flat distribution",
        7: "Galicher et al. 2016 powerlaw (a^-0.61)", #appropriate for planets
--Aout_max     {--a_max}    upper limit for the outer semi-major axis [5e6 RSun]
--Aout_min     {--a_min}    lower limit for the outer semi-major axis [0.5 RSun]
--Aout_distr   {--a_distr}  outer semi-major axis option: 
        0: "Log Uniform distribution", #default
        1: "Constant semi-major axis",
        2: "Tokovinin lognormal mu = 10^5d, sigma = 2.3",
        3: "Lognormal mu = 10^3.5d, sigma = 2.3",
        4: "Rizzuto Lognormal mu = 10^0.95 AU, sigma = 1.35",
        5: "Sana et al. 2012",
        6: "flat distribution",
        7: "Galicher et al. 2016 powerlaw (a^-0.61)", #appropriate for planets
--Ein_max      {--E_max}    upper limit for the inner eccentricity [1.]
--Ein_min      {--E_min}    lower limit for the inner eccentricity [0.]
--Ein_distr    {--E_distr}  inner eccentricity option: 
        0: "Thermal", #default
        1: "Constant eccentricity",
        2: "Sana et al. 2012 e^-0.45", #-> close binaries
        3: "Flat distribution",
        4: "Powerlaw e^0.5",   
        5: "Bowler et al. 2020 Beta distribution", #appropriate for planets                                                                          
--Eout_max     {--e_max}    upper limit for the outer eccentricity [1.]
--Eout_min     {--e_min}    lower limit for the outer eccentricity [0.]
--Eout_distr   {--e_distr}  outer eccentricity option: 
        0: "Thermal", #default
        1: "Constant eccentricity",
        2: "Sana et al. 2012 e^-0.45", #-> close binaries
        3: "Flat distribution",
        4: "Powerlaw e^0.5",
        5: "Bowler et al. 2020 Beta distribution", #appropriate for planets                                          
--i_max                     upper limit for the relative inclination [pi]
--i_min                     lower limit for the relative inclination [0]
--i_distr                   relative inclination option: 
        0: "Circular uniform distribution", #default
        1: "Constant inclination",
--Gin_max      {--G_max}    upper limit for the inner argument of pericenter [pi]
--Gin_min      {--G_min}    lower limit for the inner argument of pericenter [-pi]
--Gin_distr    {--G_distr}  inner argument of pericenter option: r
        0: "Uniform distribution", #default
        1: "Constant argument of pericenter",    
--Gout_max     {--g_max}    upper limit for the outer argument of pericenter [pi]
--Gout_min     {--g_min}    lower limit for the outer argument of pericenter [-pi]
--Gout_distr   {--g_distr}  outer argument of pericenter option: 
        0: "Uniform distribution", #default
        1: "Constant argument of pericenter",     
--Oin_max     {--O_max}    upper limit for the inner longitude of ascending node [pi]
--Oin_min     {--O_min}    lower limit for the inner longitude of ascending node [-pi]
--Oin_distr   {--O_distr}  inner longitude of ascending node option: 
        0: "Circular niform distribution", 
        1: "Constant longitude of ascending nodes", #default
        (outer longitude of ascending nodes = inner - pi)             

-T or -t                     binary end time. [13500 Myr]
-Z              {-z}         metallicity of stars  [0.02 (Solar)] 
-n                           number of triples to be simulated  [1]
-N                           number of initial triple  [0]


--SN_kick_distr   supernova kick distribution   default = 10
        options:
        0:  No kick 
        1:  Hobbs, Lorimer, Lyne & Kramer 2005, 360, 974  
        2:  Hobbs, Lorimer, Lyne & Kramer 2005, 360, 974  scaled down for bh
        3:  Arzoumanian ea 2002, 568, 289
        4:  Arzoumanian ea 2002, 568, 289 scaled down for bh
        5:  Hansen & Phinney 1997, 291, 569
        6:  Hansen & Phinney 1997, 291, 569 scaled down for bh
        7:  Paczynski 1990, 348, 485
        8:  Paczynski 1990, 348, 485 scaled down for bh
        9:  Verbunt, Igoshev & Cator, 2017, 608, 57
        10:  Verbunt, Igoshev & Cator, 2017, 608, 57 scaled down for bh #default

action items                    add these to:
--no_stop_at_merger             avoid stopping the simulation after a merger
--no_stop_at_inner_collision    avoid stopping the simulation after a collision in the inner binary
--no_stop_at_outer_collision    avoid stopping the simulation after a collision involving the outer star
--no_stop_at_disintegrated      avoid stopping after the system disintegrated into seperate systems
--no_stop_at_mass_transfer      avoid stopping the simulation at the onset of mass transfer 
--stop_at_semisecular_regime    to stop the simulation if the sytem is in the semi secular regime
--stop_at_SN                    to stop the simulation when a supernova occurs
--stop_at_CPU_time              to stop the simulation when the computational time exceeds a given value

```
## TRES-Exo for exoplanet research

In Columba et al. 2023, (A&A, 675A, 156C) we presented an extension for TRES to incorporate exoplanets.  
For exoplanet research, we have included the following processes and recommend the following settings:

1) First things first: tell TRES to include Sub-Stellar Objects by setting the boolean EXCLUDE_SSO to False in TRES.options.py. Amongst other things, this will automatically set the minimum mass of a body to 0.2 Jupiter masses (~0.0002 Solar masses) in stead of 0.0075 Solar masses. Note that our extension is only valid for giant planets, hence the minimum mass of 0.2 Jupiter masses.
2) For the dynamical stability of triples, you can which stability criterium is used through the parameter stability_limit_specification in TRES_setup.py. Options applicable to exoplanets are the simple prescription from Petrovich et al. (1), the full prescription from Petrovich et al. (2), Holman's prescription for S-type orbits (3), Holman's prescription for P-type orbits (4). 
3) energy-limited atmospheric photoevaporation of planets
4) initial planetary spin rate as 0.126 the breakup speed (Bryan+2018). This can be adjusted in the function initial_angular_frequency in triple_class.py

Examples 
1) To evolve a single system (binary star + CBP) that survives for one Hubble time (13.5 Gyr) with the following parameters:
```
- primary mass:   M1           1.04  Msun 
- secondary mass: M2           1.00  Msun
- CBP mass:       M3           0.011 Msun
- inner binary semimajor axis:  Ain       54.6  Rsun
- CBP semimajor axis:           Aout      477.8 Rsun
- inner binary eccentricity:    Ein       0.67
- CBP eccentricity:             Eout      0.15
- relative orbital inclination: i         1.9 (rad)
- simulation time:      T    1000 Myr
```

you can run TRES as:
```
python TRES.py --M1 1.04 --M2 1. --M3 0.011 --Ain 54.6 --Aout 477.8 --Ein 0.67 --Eout 0.15 -i 1.9 -T 1000  --no_stop_at_mass_transfer 
```
and use the rdc_TRES.py and rdc_TRES_csv.py to read the output and print it in readable text format. If you wish to run the complete evolution, use -T 13500 to simulate one Hubble time and you'll obtain a DWD-orbiting CBP ('Magrathea' planet).

2) To evolve a single system (binary star + CBP) where the inner binary merges as a DWD around 10 Gyr, you can use the following command:

```
python TRES.py --M1 1.33 --M2 1.06 --M3 0.0046 --Ain 26.35 --Aout 3012.9 --Ein 0.3 --Eout 0.1 -i 1.7 -T 11000  --no_stop_at_mass_transfer -f 'testRun_2.hdf'
```

## TRES with MESA

TRES can now be run with MESA as stellar code. To do so, you simply need to switch option 'USE_MESA_AS_STELLAR_CODE' in TRES_options.py to True. In this case, we also advise to switch to True 'GET_GYRATION_RADIUS_FROM_STELLAR_CODE' and 'GET_AMC_FROM_STELLAR_CODE'. These physical quantities are in this case obtained directly from the structure of the star.

If you want to change settings within MESA, this can be done in AMUSE through "particles.set_control('name_of_control', value)". For TRES, we provide the script MESA_setup.py which already sets a few MESA controls. The choice of input physics was made in order to be the closest to the SeBa defaults. More details can be found in Sciarini et al. 2025 (in prep.).

The list of MESA controls can be found in https://github.com/MESAHub/mesa/blob/r15140/star/defaults/controls.defaults.    


## Understanding the TRES output

Normally TRES adds the evolution history of individual triples in the TRES.hdf file. Every snapshot represents a moment in the evolution of the triple when something interesting happened, for example one of the star transitions from the main-sequence to the hertzsprung gap, or mass transfer starts or stops. 


## Reducing the TRES output

The python script rdc_TRES_csv.py reduce the TRES hdf output and creates a csv file. The full list of available options is [default]:
```
-f      Root of the name of the input file [TRES]
-P 	Parameter_style [1]
-S      Printing style [0] 
-F      Print all snapshots. By default only the first & last lines are printed.
```

For the parameter style you can choose between:
```
0      All parameters
1	Selected parameters [default]
2	Parameters can be selected by the user in the function create_snapshot_for_dict
```

You can also print the parameters to the terminal through the printing style option:
```
0      No printing to screen
1	Printing to screen
```



You can also select specific types of triples. For these a single extra line is added on the first occasion the requirements are met. Options are:
```
--st1          stellar type of inner binary primary star [-1]
--st2          stellar type of inner binary secondary star [-1]
--st3          stellar type of outer star [-1]
--btin         binary type of inner binary [-1]
--btout        binary type of outer binary [-1]
--trt          triple type [-1]
```
or if you prefer to specify these in string format:
```
--st1str       stellar type of inner binary primary star [all]
--st2str       stellar type of inner binary secondary star [all]
--st3str       stellar type of outer star [all]
--btinstr      binary type of inner binary [all]
--btoutstr     binary type of outer binary [all]
--trtstr       triple type [all]
```


The depreciated script rdc_TRES.py only prints to screen. Here the printing style options are:

Which parameters are printed and in which style can be adjusted to your liking in the function rdc().
Currently there are 3 options settable on the command line via -S (print_style):
```
0      Selected parameters are printed in a human readible way
1      Full - all possible parameters are printed (sys.exit after first snapshot)
2      TRES standard - selected parameters
3      TRES standard - selected parameters, csv format 
```



For option 2: 
6 lines are printed for every snapshot. The columns represent:

General information on the system: 
```
Line 1: snapshot number, triple number, time, relative_inclination, dynamical_instability, kozai_type, error_flag_secular, CPU_time
```
Orbital information (inner binary | outer binary) :
```
Line 2: 'bs:', binary type, semimajoraxis, eccentricity, argument_of_pericenter, longitude_of_ascending_node 
        | binary type, semimajoraxis, eccentricity, argument_of_pericenter, longitude_of_ascending_node 
```
Stellar information (primary | secondary | tertiary)
```
Line 3: 'st:', is_donor, stellar_type, mass, spin_angular_frequency, radius, core mass
        | is_donor, stellar_type, mass, spin_angular_frequency, radius, core mass
        | is_donor, stellar_type, mass, spin_angular_frequency, radius, core mass
```        


For option 0: 
One line is printed for every snapshot with the parameters in the same order as above (excluding the snapshot number). The units are Solar Mass, Solar radius, Myr. 

The stellar types in TRES follow the standard terminology of AMUSE:
```
0   lm_ms     deeply or fully convective low mass MS star
1   ms        Main Sequence star
2   hg        Hertzsprung Gap
3   rgb       First Giant Branch
4   cheb      Core Helium Burning
5   eagb      Early Asymptotic Giant Branch
6   tpagb     Thermally Pulsating Asymptotic Giant Branch (not used in SeBa -> labelled as 5) 
7   hems      Main Sequence Naked Helium star
8   hehg      Hertzsprung Gap Naked Helium star
9   hergb     Giant Branch Naked Helium star
10  hewd      Helium White Dwarf
11  cowd      Carbon/Oxygen White Dwarf
12  onewd     Oxygen/Neon White Dwarf
13  ns        Neutron Star
14  bh        Black Hole
15  sn        Massless Supernova
16  unknown   Unknown stellar type
17  prems     Pre-main-sequence Star
18  planet    Planet
19  bd        Brown Dwarf
```
When selecting for stellar type, you can also use the following strings for combinations of startypes:
```
agb     eagb, tpagb             [5,6]
heg     hehg, hergb             [8,9]
wd      hewd, cowd, onewd       [10,11,12]            
```
Note that stellar type 0, 6 & 8 are not used in SeBa. All main-sequence stars are stellar type 1, all AGB stars have the label stellar type 5 (both EAGB as TPAGB stars), and both helium hertzsprung gap stars as helium giants are labelled with stellar type 9.  


The binary type is a classification for a specific orbit, e.g. the inner or the outer orbit of a triple. The following options exist:
```
-1  all
0   unknown
1   merger
2   disintegrated
3   dynamical_instability
4   detached
5   contact
6   collision
7   semisecular
8   rlof
9   stable_mass_transfer
10  common_envelope
11  common_envelope_energy_balance (i.e. alpha-CE)
12  common_envelope_angular_momentum_balance (i.e. gamma-CE)
13  double_common_envelope
```

And similarly for the triple as a whole:
```
-1  all
0   hierarchical
1   dynamical_instability 
2   semisecular_regime (currently not in use)
3   error_flag_secular 
```

Do you want to rerun a system in your datafile? No need to copy all the parameters, simply run rdc_TRES.py with two extra parameters: 

```
--print_init      to print initial conditions for re-running 
-l                the line number of the first line in your hdf datafile where the system appears
                  where the stars are on the zero-age main-sequence. 
``` 
For example: ```rdc_TRES.py -f TRES.hdf --print_init -l 0```. This will return something like:
```amuse TRES.py -M 1.3 -m 0.5  -l  0.5 -A 200.0 -a 20000.0 -E 0.1 -e 0.5 -G 0.1 -g 0.5 -I 1.3962634016 ```


## TRES-development-team
For more advanced tips, see the README in the developer-folder.

## References

See the following publication: [Toonen et al 2016](https://ui.adsabs.harvard.edu/abs/2016ComAC...3....6T/abstract) for more details on TRES in general.
See [Columba et al 2023](https://ui.adsabs.harvard.edu/abs/2023A%26A...675A.156C/abstract) for more details on TRES-Exo.




