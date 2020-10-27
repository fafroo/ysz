
           

                                              

### Conductivities #####
#
# OLD_MONO_100 -> [700, 750, 800, 850] = [ ???, ???,  3.58  , 5.75]
#
#
#
#

#######################################################################################################################################
################################################# separate ############################################################################
#######################################################################################################################################                                          
                                          
############################################################################
#################### EEE ###################################################
############################################################################

################# no_e_fac ###################
# EEE EIS CV ... relatively good fit for both CV and EIS
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=100, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                         prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                          prms_values=(true, 1.0, 1.0, 1.0, 21.94051183287039, 21.553329968593776, 20.965273571151613, 0.09841737344413018, -0.091875316601428, 0.04652480321433385, 9.3e-11, 0.85, 6.478278331551995, 7.546431173856936, 0.0)
                         ,use_experiment=true);

# EEE EIS CV .. fitnes factos {1, 10} ... good fit for both ... only 2 peaks
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=100, bias=0.0, data_set="OLD_MONO_100", pyplot=1,   
                  prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                   prms_values=(true, 1.0, 1.0, 1.0, 24.383816208538917, 21.94866414517787, 20.68161992890677, 0.1380699169577035, -0.22933373067502555, -0.23743440380456354, 9.3e-11, 0.85, 2.7481343742171562, 2.930552174417333, 0.0)
                  ,use_experiment=true);                         
                         
### COOL ... without e_fac maybe the best EEE EIS CV fit                        
# EEE EIS CV ... fitness factors {1, 1} ... very good EIS fit ... 3 peaks     ######################## COOOL !!!!! ### but try pO2 = 80 !!!!                        
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=100, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                  prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                   prms_values=(true, 1.0, 1.0, 1.0, 21.96912506564684, 21.569512049276458, 20.93977987957153, 0.1241112223034025, -0.10660629726868907, -0.006580383278506099, 9.3e-11, 0.85, 5.848806563958082, 8.27404376722011, 0.0)
                  ,use_experiment=true);
                  
# EEE EIS CV ... fitness factors {0.2, 1} ... 2 peaks (+ 1 very small)      
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                  prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                   prms_values=(true, 1.0, 1.0, 1.0, 22.09624010786229, 21.66950610327097, 20.565105499265908, 0.1705724661650871, -0.1446619735003291, 0.034265126988040316, 9.3e-11, 0.85, 2.6704017916319613, 2.415004897727594, 0.0)
                  ,use_experiment=true);
                  
################## with e_fac ################# 

#!!!!!!! the following try fit with CV ... for pO2=80
# EEE EIS .... CV fitted automatically :) ... very good fit for .. with e_fac !!                  
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.CV_simulation(800, [100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                         prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                          prms_values=(true, 1.0, 1.0, 1.0, 21.9445, 21.3697, 20.8594, 0.101682, -0.178855, 0.115187, 9.3e-11, 0.85, 6.22922, 5.75616, 0.138391)
                         ,use_experiment=true);
                         
                         
                  

############################################################################
#################### LLL ###################################################
############################################################################


################# no_e_fac ##############
# 4 EIS ... not capturing the non-monotony ... even with extended search area ##### # # ## ## ## # # # 
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80, 100], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                   prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC"     ], 
                                   prms_values=(0.0, 0.0, 0.0, 23.536793805295844, 26.44428245339161, 21.308492340939917, 0.5500221642340698, 0.029112401256050482, -0.020889037864824146, 9.3e-11, 0.85, 2.0593360762011694)
                                   ,use_experiment=true);                                        

# 4 EIS ELE ... example ... but very simmilar to"separate_vacancy", LLL                                 
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80, 100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                          prms_names["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "e_fac"], 
                                          prms_values=(0.0, 1.0, 0.0, 24.4631, 21.7136, 21.5347, -0.0113085, -0.0926479, 0.179689, 9.3e-11, 0.85, 3.70135, 0.0)
                                          ,use_experiment=true);
                                          
# 4 EIS EEL ... but pretty much the same as the above LEL                                          
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80, 100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                          prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "e_fac"], 
                                          prms_values=(1.0, 1.0, 0.0, 24.1783, 21.7137, 21.5419, -0.0107187, -0.0862911, 0.185543, 9.3e-11, 0.85, 3.74939, 0.0)
                                          ,use_experiment=true);
                    
                    
################# with_e_fac ##############                    
                 
# EEL EIS .... CV hooked ... sensitive on pO2 ... 3 peaks... with e_fac
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40,60,80, 100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                         prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                          prms_values=(1.0, 1.0, 0.0, 21.9258, 21.2852, 21.7927, 0.0277078, -0.280821, 0.362658, 9.3e-11, 0.85, 7.54141, 1.65368, 0.296044)
                         ,use_experiment=true);
                         
# LEE EIS .... CV hooked
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40,60,80, 100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(0.0, 1.0, 1.0, 24.253505704719995, 27.895308798638613, 20.689058014097917, 0.5435406241193925, 0.2124972805439068, -0.025572361570387683, 9.3e-11, 0.85, 2.1001141086263937, 0.8016033006905938, 0.004387336514558263)
                                ,use_experiment=true);

# ELE EIS .... CV hooked ... the worse sensitivity to pO2 ... and the worse shape of Nyquisty
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40,60,80, 100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(1.0, 0.0, 1.0, 22.004601154444064, 21.969771021412477, 20.64499552886603, -0.5553757496916619, -0.04999999999933978, 0.23448570623706932, 9.3e-11, 0.85, 5.925253670629901, 9.734714051365351, 0.7504653628307916)
                                ,use_experiment=true);                 
                 
#### LLL CV .... husty!!! ... trefa .... ale EIS je tragicke
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.CV_simulation(800, [100], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(true, 0.0, 0.0, 0.0, 24.2177, 21.6481, 22.6917, 0.698797, 0.252055, -0.122604, 9.3e-11, 0.85, 1.46538, 4.51906, 0.285749)
                                ,use_experiment=true);                    
##########
## COOL ##                                 
#### LLL EIS CV .... husty!!! ... obe nejsou spatni                                   
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=[40, 60, 80], bias=0.0, data_set="OLD_MONO_100", pyplot=1,   
                                              prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                               prms_values=(true, 0.0, 0.0, 0.0, 22.576839213753235, 22.011902293420093, 27.69994335225577, 0.04823267858003028, -0.2710822686942347, 0.5656693158734294, 9.3e-11, 0.85, 0.21671402944207255, 9.144064551170423, 0.3033398196955781)
                                              ,use_experiment=true);                    
                    
### LLL EIS CV .... BEZ e_fac   # --->> nic moc :)
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                           prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                           prms_values=(true, 0.0, 0.0, 0.0, 24.069480130382026, 27.84830853588366, 22.089454852605634, 0.7999985765425789, -0.5103390865926676, -0.022909285494503932, 9.3e-11, 0.85, 0.9567488020824775, 2.851374822442069, 0.0)
                           ,use_experiment=true);                    
                    
                    
s

                                          
                                          
                                          
                                          
### HAND FIT ### :) LLL EIS ... quite good EIS... bad CV
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 
                                                     22.277, 22.312, 22.3, 
                                                     0.048, -0.271, 0.566, 
                                                     1.05e-11, 0.35, 24., 2.544, 0.303)
                                       ,use_experiment=true);
                                       
#                                                    #################################################
### corrected HAND FIT ... LLL EIS ... WOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOW WOWOWOWWOOW WO !!!!!! BRUTALISK !!!!! EIS cooool !!!!
#                                                    #################################################
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 22.325513444696398, 22.100276559559344, 22.382691084143286, 0.04243114428504254, -0.22279524777844634, 0.5430896356175956, 1.05e-11, 0.35, 34.397339021979015, 4.503994521917953, 0.338333333916094)
                                       ,use_experiment=true);
                                       
######### corrected HAND FIT for more pO2 ... LLL EIS 
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 22.551514573009637, 22.221241593154982, 22.514061660138903, 0.03948393867277375, -0.03391657645569906, 0.6777921770403312, 1.05e-11, 0.35, 47.67570795526793, 3.363950713557616, 0.6107672677397736)
                                       ,use_experiment=true);
#### .... continuation from previos... to pO2 60
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 22.519825220050315, 22.101114333751948, 22.537917125481112, 0.046779285103424426, -0.036495661417036834, 0.6645252416469626, 1.06e-11, 0.35, 66.81244689470016, 2.899812097038787, 0.7357178655306754)
                                       ,use_experiment=true);
                           
##### ....... and ... then pO2 40, 60, 80,  again                           
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 22.49042860963425, 22.18758070122782, 22.535374543362902, 0.04626751671867485, -0.02797383209333109, 0.6910672027823139, 1.06e-11, 0.35, 51.10551528037439, 3.2209400174059213, 0.6775179360491119)
                                       ,use_experiment=true);
                    

#!#!!##!!#!#!#!#!#!#!#!#!#!#!#!#!!##!!#!#!#!#!#!#!##!#!!##!!#!#!#!#!#!#!#!#!#!#!#!#!!##!!#!#!#!#!#!#!#
### Temperature study !!#!#!#!#!!##!!#!#!#!#!#!#!#!#!#!#!#!#!!##!!#!#!#!#!#!#!#

#### only EIS fit
# 700 ... 2 EIS LLL ...     EIS very good even for pO2 20 !!! ... but CV bad
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(700, [20, 40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(1, 0.0, 0.0, 0.0, 22.0248, 21.6486, 21.5078, -0.12619, 0.0824053, -0.129584, 1.02, 0.752653, 9.66799, 7.46362, 0.020354)
                                ,use_experiment=true); 
                                
# 750 ... 2 EIS LLL ...    EIS great even for pO2 20 !!! ... but CV bad
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(750, [20, 40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(1, 0.0, 0.0, 0.0, 22.3038, 22.0987, 22.0554, -0.209048, -0.0667363, 0.239006, 2.07, 0.38203, 27.384, 11.7518, 0.0148614)
                                ,use_experiment=true);

# 800 ... 2 EIS LLL ... 
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 22.551514573009637, 22.221241593154982, 22.514061660138903, 0.03948393867277375, -0.03391657645569906, 0.6777921770403312, 1.05e-11, 0.35, 47.67570795526793, 3.363950713557616, 0.6107672677397736)
                                       ,use_experiment=true); 
    #alternative  
    new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(true, 0.0, 0.0, 0.0, 22.56783400557656, 22.22237966317968, 22.50693968907697, 0.05520337660564685, -0.033209183842112, 0.6770403513536281, 3.62, 0.35, 45.08895790653341, 3.269908550672186, 0.6182027814202923)
                                ,use_experiment=true);  
    
    # alternative #2
    new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(true, 0.0, 0.0, 0.0, 23.0467, 22.5346, 22.4665, -0.05, 0.0525283, 0.697041, 9.3e-11, 0.849408, 38.7172, 7.74472, 0.651215)
                                ,use_experiment=true); 
                                
    #alternative  #3 .. new one
    new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                 prms_values=(true, 0.0, 0.0, 0.0, 23.671342687700914, 22.75009386895551, 22.36278147296789, -0.07877895777799904, 0.10251592543765028, 0.633939269456473, 3.62, 0.8503249682155306, 29.210203301013003, 12.194786069416887, 0.416244601465356)                                
                                ,use_experiment=true);  
                                
                                
                                
# 850 ... 2 EIS LLL ... 
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(850, [40, 60], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 25.474103553624037, 22.925019737001634, 21.85804785753094, -0.10831702948372236, 0.38984727556707693, -0.23477140773747657, 5.85, 0.65, 3.4114628080433516, 20.67717563382255, 0.0)
                                       ,use_experiment=true);
                             
           
#### EIS-CV fits


## 700 ... 2 EIS-CV LLL  ... ff1-100
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(700, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,   
                                              prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                               prms_values=(1, 0.0, 0.0, 0.0, 
                                                            21.9222, 21.2101, 21.809,
                                                            0.0558328, 0.0480154, -0.0438488, 
                                                            1.02, 0.765282, 37.9658, 5.55807, 0.270473)
                                              ,use_experiment=true);
           
## 750 ... 3 EIS-CV LLL ... ff ???
# first attempt with ff1-100 --> EIS too bad
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(750, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,
                                              prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                               prms_values=(1, 0.0, 0.0, 0.0, 22.2922, 21.4767, 26.0711, 0.0326867, -0.0220428, 0.0672773, 2.07, 0.83353, 34.7543, 7.64352, 0.251502)
                                              ,use_experiment=true);
 
# second attempt with ff5-100 --> much better, 3 peaks
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=750, pO2=[40, 60], bias=0.0, data_set="OLD_MONO_100", pyplot=1,   
                                                     prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                                      prms_values=(1, 0.0, 0.0, 0.0, 
                                                           22.178, 21.508, 22.1441, 
                                                           0.00928035, 0.210568, -0.104948, 
                                                           2.07, 0.475152, 41.6055, 9.92313, 0.339722)
                                                     ,use_experiment=true);

# third attempt with ff10-100 -->  ok .. :) we can take it
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=750, pO2=[40, 60], bias=0.0, data_set="OLD_MONO_100", pyplot=1,   
                                                     prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "conductivity", "nu", "OC", "ms_par", "e_fac"],
                                                      prms_values=(1, 0.0, 0.0, 0.0, 
                                                          22.1636, 21.4919, 22.0697,  
                                                          -0.00599031, 0.18108, -0.05,  
                                                          2.07, 0.532246, 51.0891, 9.42743, 0.362266)
                                                     ,use_experiment=true);                                                     
                                                     

## 800 ... 2 EIS_cv LLL ... ff ???
# some previous attempt ... unknown ff ... only 2 peaks
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=[40, 60, 80], bias=0.0, data_set="OLD_MONO_100", pyplot=1,   
                                              prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                               prms_values=(true, 0.0, 0.0, 0.0, 
                                               22.576839213753235, 22.011902293420093, 27.69994335225577, 
                                               0.04823267858003028, -0.2710822686942347, 0.5656693158734294, 
                                               9.3e-11, 0.85, 0.21671402944207255, 9.144064551170423, 0.3033398196955781)
                                              ,use_experiment=true);
                                              
                    
# sedonc
(1, 0.0, 0.0, 0.0, 22.4555, 21.6669, 22.7248, -0.0285608, 0.195182, 0.102502, 3.58, 0.0502756, 28.2008, 17.8957, 0.4)




## 850 ... 2 EIS_cv LLL ... ff ???
           
           
           
           
           
           
           
           
           
           
# fit pres bety a S ... nic moc stejne! :( 
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=700, pO2=[40, 60, 80], bias=0.0, data_set="OLD_MONO_100", pyplot=1, 
                                prms_names=["separate_vacancy",
              "A.beta", "R.beta", "O.beta",
              "A.S", "R.S", "O.S",
              "A.exp", "R.exp", "O.exp", 
              "A.r", "R.r", "O.r",              
              "A.DG", "R.DG", "O.DG",     
              "conductivity", "nu",      "OC", "ms_par", "e_fac"  ],
                                 prms_values=(1, 0.5468986658625657, 0.47998951044726024, 0.2779139613165412, -0.27769563747657716, -0.09593075662029801, -0.22320858430182988, 0.0, 0.0, 0.0, 21.9464, 21.5159, 21.4928, 0.0032943, 0.0277383, -0.0838499, 1.02, 0.811907, 10.472, 6.00553, 0.0718344)
                                ,use_experiment=true);             
           
           
           
           
           
#######################################################################################################################################
################################################# shared ##############################################################################
#######################################################################################################################################

############################################################################
#################### EEE ###################################################
############################################################################

                                          
##### MAYBEE...... another hand-fit 
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(false, 1.0, 1.0, 1.0, 
                                                            21.9, 21.3, 21.0, 
                                                            0.15, 0.0, -0.40967, 
                                                            1.06e-11, 0.35, 1.0, 20., 0.0)
                                                ,use_experiment=true);
############
### COOL ###
### EEE 1 EIS ... sha ... hustYYY!! ...                                                 
ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(false, 1.0, 1.0, 1.0, 21.78592822930247, 21.347238180487707, 20.845011931720087, 0.11415904394451544, -0.04710867464932109, -0.4221231346344268, 1.06e-11, 0.35, 1.252221527673722, 37.34152581968751, 0.0)
                                                ,use_experiment=true);                                                
                                                
                           
### EEE 1 EIS-CV ... sha ... fitness = {1 ,1}   --->> great EIS, medium CV
ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(false, 1.0, 1.0, 1.0, 21.77097654529923, 21.339094620305584, 20.83825101342013, 0.11612661311085312, -0.01621306293886935, -0.44405718790607684, 1.06e-11, 0.35, 0.6404908464083954, 37.19354903267066, 0.010059997045145626) 
                                                ,use_experiment=true);
                           
### EEE 1 EIS-CV ... sha ... fitness = {0.5 ,1} --->> very good EIS, not good CV                          
ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(false, 1.0, 1.0, 1.0, 21.768338176697593, 21.339959521197464, 20.77053224186594, 0.14211354657568556, -0.06522415408171552, -0.4099634528362534, 1.06e-11, 0.35, 0.9838424541352313, 43.614438009170286, 0.04296470330868719)
                                                ,use_experiment=true);                           
                           
### EEE 1 EIS-CV ... sha ... fitness = {0.1 ,1} --->> EIS 2 peaks, good CV                         
ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(false, 1.0, 1.0, 1.0, 27.204738430905046, 21.042465727359954, 20.903492276526453, -0.7650977954784123, -0.16483434129293068, -0.5474897673095572, 1.06e-11, 0.35, 7.553853360574966, 20.131456245497063, 0.12330883448476397)
                                                ,use_experiment=true);                            
                           
### EEE EIS CV.... shared vac ... not bad, not the best ... 2 peaks
ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                          prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC",  "ms_par","e_fac"     ], 
                                          prms_values=(false, 1.0, 1.0, 1.0, 21.6797, 21.0495, 20.7467, 0.0212226, -0.423638, 0.798779, 9.3e-11, 0.85, 6.75588, 9.95628, 0.305956)
                                          ,use_experiment=true);                                        
                          
### EEE 1 CV ... sha ...   --->> EIS 1 peak, better CV                         
ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                                prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(false, 1.0, 1.0, 1.0, 27.280316820837967, 20.744018794803896, 25.703864047655617, -0.10168151776708975, -0.10239784764049484, -0.7893449634234231, 1.06e-11, 0.35, 0.0029974884900697343, 8.174364689405985, 0.11692539726030851)                           
                                                ,use_experiment=true);                                           
                                          

############################################################################
#################### LLL ###################################################
############################################################################


#### LLL 1 EIS .... shared vac   -> CV bad , EIS 2 peaks

#### ...... LLL shared ... 1 EIS CV ... fitness = {10, 10}
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                              prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                               prms_values=(0, 0.0, 0.0, 0.0, 22.6529, 21.2666, 22.653, 0.201058, 0.0140769, -0.223763, 1.06e-11, 0.35, 1.33437, 25.0502, 0.644818)
                                              ,use_experiment=true);
                                              
#### LLL shared ... 1 EIS CV .. fitness = {1, 10}
new_GAS_LoMA_shared = ysz_fitting.simple_run(simulations=["CV", "EIS"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                                     prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                                      prms_values=(false, 0.0, 0.0, 0.0, 23.660999804374093, 21.16644727635489, 23.246774488288892, 0.5590363181503285, 0.3057035105534315, -0.46167785857357707, 1.06e-11, 0.35, 57.20049929904167, 10.949218885809483, 0.40238189628810334)
                                                     ,use_experiment=true);                                              

#### LLL shared ... 1 EIS CV .. fitness = {1, 5}                                                     
ysz_fitting.simple_run(simulations=["EIS", "CV"], TC=800, pO2=40, bias=0.0, data_set="OLD_MONO_100", pyplot=1,  
                                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                       prms_values=(0, 0.0, 0.0, 0.0, 22.3718, 22.159, 26.9917, 0.00119476, -0.217313, 0.74843, 1.06e-11, 0.35, 12.344, 11.4685, 0.383988)
                                                       
                                                       ,use_experiment=true);                                                     

#### ...... LLL shared !!! 2 EIS CV ...not bad :)   >>>>>>> almost THE SAME result as for the separate vac
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                              prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                               prms_values=(false, 0.0, 0.0, 0.0, 22.5615, 22.113, 27.0515, 0.0456777, -0.33161, 0.684794, 9.3e-11, 0.85, 75.0419, 9.5712, 0.311755)
                                              ,use_experiment=true);



                                                                         

                           
                           
                           
                           
                           
                           
                           
                           
        
                           
                           
                            

                           

                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
###########################################################                           
#####################  PAR_STUDY  #########################                           
# separate VS shared vacancy ... 3 vs 2 peaks
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                           prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                           prms_values=(true, 1.0, 1.0, 1.0, 
                                       21.6106, 21.898, 21.4851, 
                                       0.497806, collect(-0.3 : -0.02 : -0.3), 0.229967, 
                                       9.3e-11, 0.85, collect(1.0 : 1 : 14), 3.5, 0.0)
                           ,use_experiment=false);

                           
# e_fac descussion ... without inductance
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                           prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                           prms_values=(0.0, true, 0.0, 0.0, 0.0, 
                                       23.6106, 22.598, 23.4851, 
                                       -0.497806, collect(-0.3 : -0.02 : -0.3), 0.229967, 
                                       9.3e-11, 0.85, 1.0, 1.0, collect(0.0 : 0.1 : 0.5))
                           ,use_experiment=false);
                           
                           
                           
                           
                           
                           
# OC study ... 3 clear peaks  LLL   !!!!                       
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                           prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                           prms_values=(0.0, [true, false], 0.0, 0.0, 0.0, 
                                       22.2106, 22.098, 21.5851, 
                                       -0.1497806, 0.0, 0.029967, 
                                       9.3e-11, 0.85, 10.0, 10., collect(0.0 : 0.1 : 0.0))
                           ,use_experiment=false); 
# separate_vacancy study ... 3 clear peaks ... LLL !!!
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                  prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                  prms_values=(0.0, [true, false], 0.0, 0.0, 0.0, 
                                              22.2106, 21.798, 21.1851, 
                                              -0.1497806, -0.0, -0.19967, 
                                              9.3e-11, 0.85, 10.0, 10., collect(0.0 : 0.1 : 0.0))
                                  ,use_experiment=false);
                           
# OC study ... 3 clear peaks ... EEE !!!!                          
                           ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                  prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                  prms_values=(0.0, [true, false], 1.0, 1.0, 1.0, 
                                              22.1106, 22.098, 21.5851, 
                                              0.1497806, -0.0, 0.029967, 
                                              9.3e-11, 0.5, 1000.0, 100., collect(0.0 : 0.1 : 0.0))
                                  ,use_experiment=false);
                           
                           
########################################################################################################                           
#  .... .. .. 4 peaks clear!... EEE !!! bulk peak!
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                                       prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                       prms_values=(0.0, [true], 1.0, 1.0, 1.0, 
                                                                   21.2, 21., 21.0, 
                                                                   0., 0., 0., 
                                                                   [0.05]*1.0e-11 , [0.9], 30., 10., collect(0.0 : 0.1 : 0.0))
                                                       ,use_experiment=false);
########################################################################################################
                           
                           
                           
# e_fac study .... 3 clear peaks ... down loping !!!! :D
# the same looping is appearing using EXP kinetics. 
#### when A or O is fast, middle peak starts to be the boundary one.... and makes a "inductance tail" 
# --- dva nebo tri obrazky
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                  prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                  prms_values=(0.0, true, 0.0, 0.0, 0.0, 
                                              22.2106, 22.098, 21.5851, 
                                              -0.1497806, 0.0, 0.029967, 
                                              9.3e-11, 0.85, 10.0, 10., collect(0.0 : 0.2 : 0.5))
                                  ,use_experiment=false);


                                  
### DD vs nu  --->> DD changes only R_ohm ... nu changes DRT
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                                prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(0.0, [true], 1.0, 1.0, 1.0, 
                                                            22.0, 22., 21.5, 
                                                            0.15, 0.0, -0.629967, 
                                                            [0.1, 1, 10, 100].*1.0e-11, collect(0.1 : 666.1 : 0.8), 1000.0, 100., collect(0.0 : 0.1 : 0.0))
                                                ,use_experiment=false); 
                                  
                                  
###########################################################################################
##########################        USED IN ARTICLE      ####################################
###########################################################################################
###########################################################################################
###########################################################################################

# par_stud_e_fac_LF
new_GAS_LoMA_shared = ysz_fitting.simple_run(
                                                       ysz_fitting.EIS_simulation(800, [80], 0.0, data_set="OLD_MONO_100", plot_legend=false), 
                                                       #ysz_fitting.CV_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), 
                                                       
                                                           pyplot=1,  
                                                            prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                                             prms_values=(0.0, 0.0, 0.0, 22.5, collect(22.9 : 0.1 : 22.9), 27.69994335225577, 0.04823267858003028, -0.2710822686942347, 0.5656693158734294, 9.3e-11, 0.85, 0.21671402944207255, 9.144064551170423, collect(0.0 : 0.05 : 0.2))
                                                            ,use_experiment=false);
                                                            
# par_stud_e_fac_HF
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                          prms_names=["A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"], 
                                          prms_values=(1.0, 1.0, 1.0, 21.6106, 27.898, 20.4851, -0.497806, -0.234529, 0.229967, 9.3e-11, 0.85, 3.04504, 0.05, [0.0, 0.1, 0.2, 0.3, 0.4])
                                          ,use_experiment=false);

                                          
# par_study_sha_sep_2_peaks >>> LaTeX tab. parameters: id 0.1
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                         prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                         prms_values=(0.0, [true, false], 1.0, 1.0, 1.0, 
                                                     22.0, 22., 21.5, 
                                                     0., 0.0, 0., 
                                                     9.3e-11, 0.5, 300., 100., collect(0.0 : 0.1 : 0.0))
                                         ,use_experiment=false);

# par_study_sha_sep_3_peaks >>> LaTeX tab. parameters: id 0.2
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                         prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                         prms_values=(0.0, [true, false], 1.0, 1.0, 1.0, 
                                                     22.0, 22., 21.5, 
                                                     0.15, 0.0, -0.629967, 
                                                     9.3e-11, 0.5, 1000.0, 100., collect(0.0 : 0.1 : 0.0))
                                         ,use_experiment=false);

# VERY interesting !!!
# separare_vacancy interchanged the order of reactions <-> peaks bijections >>> LaTeX tab.: id 0.2
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                                prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(0.0, [true, false], 1.0, 1.0, 1.0, 
                                                            22.0, [22., 29.], 21.5, 
                                                            0.15, 0.0, -0.629967, 
                                                            9.3e-11, 0.5, 1000.0, 100., collect(0.0 : 0.1 : 0.0))
                                                ,use_experiment=false);

# par_study_DGR >>> LaTeX tab. parameters: id 0.1                                         
ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40], 0.0, data_set="OLD_MONO_100"), pyplot=1,  
                                                prms_names=["L", "separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG",      "DD", "nu", "OC", "ms_par", "e_fac"],              
                                                prms_values=(0.0, [true], 1.0, 1.0, 1.0, 
                                                            22.0, 22., 21.5, 
                                                            0., collect(-0.3 : 0.1 : 0.1), 0., 
                                                            9.3e-11, 0.5, 300., 100., collect(0.0 : 0.1 : 0.0))
                                                ,use_experiment=false);
                                                
# HOT candidate !!! corrected_hand_guess -> expanded to 3 pO2 ..... CV rubbish
new_GAS_LoMA_shared = ysz_fitting.simple_run(ysz_fitting.EIS_simulation(800, [40, 60, 80], 0.0, data_set="OLD_MONO_100", plot_legend=false), pyplot=1,  
                                       prms_names=["separate_vacancy", "A.exp", "R.exp", "O.exp", "A.r", "R.r", "O.r", "A.DG", "R.DG", "O.DG", "DD", "nu", "OC", "ms_par", "e_fac"],
                                        prms_values=(true, 0.0, 0.0, 0.0, 22.551514573009637, 22.221241593154982, 22.514061660138903, 0.03948393867277375, -0.03391657645569906, 0.6777921770403312, 1.05e-11, 0.35, 47.67570795526793, 3.363950713557616, 0.6107672677397736)
                                       ,use_experiment=true);
