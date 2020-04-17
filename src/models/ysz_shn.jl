module ysz_shn

#############################################
# WHILE CREATING NEW MODEL                  #
# DO NOT FORGET TO CHECK                    #
#
# [x] boundary conditions
# [x] initial conditions
# [x] equilibrium phi
# [x] meas_tran & meas_stdy
# [x] output species, names
# [ ] function set_parameters
# [ ] changed module name :)

using Printf
using VoronoiFVM
using PyPlot
using DataFrames
using CSV
using LeastSquaresOptim
using ForwardDiff

const bulk_species = (iphi, iy) = (1, 2)
const surface_species = (iyAs, iyOs, iyOmins) = (3, 4, 5)
const surface_names = ("yAs", "yOs", "yOmins")

const index_driving_species = iphi

mutable struct YSZParameters <: VoronoiFVM.AbstractData

    # adsorption from YSZ
    rA::Float64   # surface adsorption coefficient [ m^-2 s^-1 ]
    DGA::Float64 # difference of gibbs free energy of adsorption  [ J ]
    betaA::Float64 # symmetry of the adsorption    
    SA::Float64 # stechiometry compensatoin of the adsorption
    kappaA::Float64 # 
    gammaA::Array{Int32,1}

    # electroreaction
    rR::Float64 # exhange current density [m^-2 s^-1]
    DGR::Float64 # difference of gibbs free energy of electrochemical reaction [ J ]
    betaR::Float64 # symmetry of the electroreaction
    SR::Float64 # stechiometry compensatoin of the electroreaction
    kappaR::Float64 # 
    gammaR::Array{Int32,1}
    
    # adsorption from gas
    rO::Float64 
    DGO::Float64
    betaO::Float64
    SO::Float64
    kappaO::Float64 # bool deciding if EXP should be used instead of LoMA
    gammaO::Array{Int32,1}
   
    # (B) electroreaction
    rB::Float64 # exhange current density [m^-2 s^-1]
    DGB::Float64 # difference of gibbs free energy of electrochemical reaction [ J ]
    betaB::Float64 # symmetry of the electroreaction
    SB::Float64 # stechiometry compensatoin of the electroreaction
    kappaB::Float64 # bool deciding if EXP should be used instead of LoMA
    gammaB::Array{Int32,1}
 
    # (C) electroreaction
    rC::Float64 # exhange current density [m^-2 s^-1]
    DGC::Float64 # difference of gibbs free energy of electrochemical reaction [ J ]
    betaC::Float64 # symmetry of the electroreaction
    SC::Float64 # stechiometry compensatoin of the electroreaction
    kappaC::Float64 # bool deciding if EXP should be used instead of LoMA
    gammaC::Array{Int32,1}
 
    stoichiometric_matrix::Array{Int32,2}

    # inductance
    L::Float64  
    
    # oxygen adsorption sites coverage w.r.t. one surface YSZ cell
    OC::Float64
    
    
    
    # fixed
    DD::Float64   # diffusion coefficient [m^2/s]
    pO2::Float64 # O2 partial pressure [bar]\z
    T::Float64      # Temperature [K]
    nu::Float64    # ratio of immobile ions, \nu [1]
    nus::Float64    # ratio of immobile ions on surface, \nu_s [1]
    numax::Float64  # max value of \nu 
    nusmax::Float64  # max value of  \nu_s
    x_frac::Float64 # Y2O3 mol mixing, x [%] 
    chi::Float64    # dielectric parameter [1]
    m_par::Float64
    ms_par::Float64

    # known
    vL::Float64     # volume of one FCC cell, v_L [m^3]
    areaL::Float64 # area of one FCC cell, a_L [m^2]


    e0::Float64
    eps0::Float64
    kB::Float64  
    N_A::Float64 
    zA::Float64  
    mO::Float64  
    mZr::Float64 
    mY::Float64
    zL::Float64   # average charge number [1]
    yB::Float64   # electroneutral value [1]
    
    phi_eq::Float64 # equilibrium voltage [V]
    y0_eq::Float64
    yAs_eq::Float64
    yOs_eq::Float64
    yOmins_eq::Float64
    
    ML::Float64   # averaged molar mass [kg]
    #
    ARO_mode::Bool # switches off B & C reactions, assumes simpler surface lattice sites
    #
    YSZParameters()= YSZParameters( new())
end

function YSZParameters(this)
    
    this.e0   = 1.602176565e-19  #  [C]
    this.eps0 = 8.85418781762e-12 #  [As/(Vm)]
    this.kB   = 1.3806488e-23  #  [J/K]
    this.N_A  = 6.02214129e23  #  [#/mol]
    this.mO  = 16/1000/this.N_A  #[kg/#]
    this.mZr = 91.22/1000/this.N_A #  [kg/#]
    this.mY  = 88.91/1000/this.N_A #  [kg/#]

    # oxide adsorption from YSZ
    this.rA= 10.0^21
    this.DGA= 0.0905748 * this.e0 # this.e0 = eV
    this.betaA = 0.5
    this.SA= 10^0.0
    this.kappaA= 0.0
    this.gammaA= [0, 1, -1, 0, 0, 0]

    # electron-transfer reaction
    this.rR= 10.0^22
    this.DGR= -0.008014 * this.e0
    this.betaR= 0.5
    this.SR= 10^0.0
    this.kappaR = 0.0
    this.gammaR= [0, 0, 1, -1, 0, 0]
    
    # (O) oxygen adsorption from gas
    this.rO= 10.0^20
    this.DGO= 0.0905748 * this.e0 # this.e0 = eV
    this.betaO = 0.5
    this.SO= 10^0.0
    this.kappaO=0.0
    this.gammaO= [0, 0, 0, 2, 0, -1]

    # (B) electron-transfer reaction
    this.rB= 10.0^20
    this.DGB= 0.001 * this.e0
    this.betaB= 0.5
    this.SB= 10^0.0
    this.kappaB = 0.0
    this.gammaB= [0, 0, 0, -1, 1, 0]

    # (C) electron-transfer reaction
    this.rC= 10.0^21
    this.DGC= (this.DGR - this.DGB)
    this.betaC= 0.5
    this.SC= 10^0.0
    this.kappaC = 0.0
    this.gammaC= [0, 0, 1, 0, -1, 0]
    
    this.stoichiometric_matrix = hcat(this.gammaA, this.gammaR, this.gammaO, this.gammaB, this.gammaC)
 
    this.L=2.3560245927364395e-6
    
    # oxygen adsorption sites coverage w.r.t. one surface YSZ cell
    this.OC = 1/4
    
    #this.DD=1.5658146540360312e-11  # [m / s^2]fitted to conductivity 0.063 S/cm ... TODO reference
    #this.DD=8.5658146540360312e-10  # random value  <<<< GOOOD hand-guess
    #this.DD=9.5658146540360312e-10  # some value  <<<< nearly the BEST hand-guess
    this.DD=4.35e-13  # testing value
    this.pO2=1.0                   # O2 atmosphere 
    this.T=1073                     
    this.nu=0.85                     # assumption
    this.nus=0.21                    # assumption
    this.x_frac=0.13                # 13% YSZ
    this.chi=27.e0                  # from relative permitivity e_r = 6 = (1 + \chi) ... TODO reference
    this.m_par =  2                 
    this.ms_par = 0.05        
    this.numax = (2+this.x_frac)/this.m_par/(1+this.x_frac)
    this.nusmax = (2+this.x_frac)/this.ms_par/(1+this.x_frac)

    this.vL=3.35e-29
    this.zA  = -2;

    this.ARO_mode = false # true# 
    #this.zL  = 4*(1-this.x_frac)/(1+this.x_frac) + 3*2*this.x_frac/(1+this.x_frac) - 2*this.m_par*this.nu
    #this.yB  = -this.zL/(this.zA*this.m_par*(1-this.nu))
    #this.ML  = (1-this.x_frac)/(1+this.x_frac)*this.mZr + 2*this.x_frac/(1+this.x_frac)*this.mY + this.m_par*this.nu*this.mO
    # this.zL=1.8182
    # this.yB=0.9
    # this.ML=1.77e-25   
    
    #phi_eq=0
    #y0_eq=0
    #yAs_eq=0
    #yOs_eq=0
    
    return this
end

# conversions for the equilibrium case
function y0_to_phi(this::YSZParameters, y0)
    yB = this.yB
    return - (this.kB * this.T / (this.zA * this.e0)) * log(y0/(1-y0) * (1-yB)/yB )
end

function y0_activity_to_phi(this::YSZParameters, a_y0)
  yB = this.yB
  return - (this.kB * this.T / (this.zA * this.e0)) * log(a_y0 * (1-yB)/yB )
end

function phi_to_y0(this::YSZParameters, phi)
    yB = this.yB
    X  = yB/(1-yB)*exp.(this.zA*this.e0/this.kB/this.T* (-phi))
    return X./(1.0.+X)
end

function equilibrium_boundary_conditions(this::YSZParameters)
    if this.ARO_mode
        a_yOs = this.pO2^(1/2.0)*exp(-this.DGO/(2 * this.kB * this.T))
        a_yAs = a_yOs*exp(-this.DGR/(this.kB * this.T))
        a_y0 = a_yAs*exp(this.DGA/(this.kB * this.T))
        return y0_activity_to_phi(this, a_y0), a_y0/(1 + a_y0), a_yAs/(1 + a_yAs), a_yOs/(1 + a_yOs), 0.2
    else   # !!!!!!!!!!!!!!!!!!!!!!!!!!
        a_yOs = this.pO2^(0.5)*exp(-this.DGO/(2.0 * this.kB * this.T))
        a_yAs = a_yOs*exp(-this.DGR/(this.kB * this.T))
        a_y0 = a_yAs*exp(-this.DGA/(this.kB * this.T))

        a_yOmins = a_yOs*exp(this.DGB/(this.kB * this.T))
        ysV = 1/(1 + a_yAs + a_yOs + a_yOmins)

        return y0_activity_to_phi(this, a_y0), a_y0/(1 + a_y0), a_yAs*ysV, a_yOs*ysV, a_yOmins*ysV
    end
end

# boundary conditions
function set_typical_boundary_conditions!(sys, parameters::YSZParameters)
    sys.boundary_values[iphi,1]=parameters.phi_eq
    sys.boundary_values[iphi,2]=0.0e-3
    #
    sys.boundary_factors[iphi,1]=VoronoiFVM.Dirichlet
    sys.boundary_factors[iphi,2]=VoronoiFVM.Dirichlet
    #
    sys.boundary_values[iy,2]=parameters.yB
    sys.boundary_factors[iy,2]=VoronoiFVM.Dirichlet
end

# initial conditions
function get_typical_initial_conditions(sys, parameters::YSZParameters)   
#function get_equilibrium_state(sys, parameters::YSZParameters)   
    inival=unknowns(sys)
    inival.=0.0
    grid = sys.grid
    #
    treshold_for_linear_function = 0.6*10^(-9)
    for inode=1:size(inival,2)
        x = nodecoord(grid, inode)[1]
        if x < treshold_for_linear_function
          inival[iphi, inode] = (((treshold_for_linear_function - x)*parameters.phi_eq)/
                                  treshold_for_linear_function)                                  
        else  
          inival[iphi,inode]=-0.0
        end
        if x < treshold_for_linear_function
          inival[iy, inode] =   (((treshold_for_linear_function - x)*parameters.y0_eq + x*parameters.yB)/
                                 treshold_for_linear_function)
        else
          inival[iy,inode]= parameters.yB
        end
    end
    inival[iyAs,1] = parameters.yAs_eq
    inival[iyOs,1] = parameters.yOs_eq
    inival[iyOmins,1] = parameters.yOmins_eq
    return inival
end


function YSZParameters_update!(this::YSZParameters)
    this.areaL=(this.vL)^0.6666
    this.numax = (2+this.x_frac)/this.m_par/(1+this.x_frac)
    this.nusmax = (2+this.x_frac)/this.ms_par/(1+this.x_frac)   
    
    this.zL  = 4*(1-this.x_frac)/(1+this.x_frac) + 3*2*this.x_frac/(1+this.x_frac) - 2*this.m_par*this.nu
    this.yB  = -this.zL/(this.zA*this.m_par*(1-this.nu))
    this.ML  = (1-this.x_frac)/(1+this.x_frac)*this.mZr + 2*this.x_frac/(1+this.x_frac)*this.mY + this.m_par*this.nu*this.mO
    
    this.phi_eq, this.y0_eq, this.yAs_eq, this.yOs_eq, this.yOmins_eq = equilibrium_boundary_conditions(this)
    return this
end

function printfields(this)
    for name in fieldnames(typeof(this))
        @printf("%8s = ",name)
        println(getfield(this,name))
    end
end

function set_parameters!(this::YSZParameters, prms_values, prms_names)
  found = false
  for (i,name_in) in enumerate(prms_names)
    found = false
    for name in fieldnames(typeof(this))
      if name==Symbol(name_in)
        if name_in in ["rA", "rR", "rO", "SA", "SR", "SO"]
          setfield!(this, name, Float64(10.0^prms_values[i]))
        elseif name_in in ["DGA", "DGR", "DGO"]
          setfield!(this, name, Float64(prms_values[i]*this.e0))   #  [DGA] = eV
        else 
          setfield!(this, name, Float64(prms_values[i]))
        end
        found = true
        break
      end
    end
    if !(found)
      println("ERROR: set_parameters: parameter name \"$(name_in)\" not found!")
      throw(Exception)
    end
  end
end





# time derivatives
function storage!(f,u, node, this::YSZParameters)
    f[iphi]=0
    f[iy]=this.mO*this.m_par*(1.0-this.nu)*u[iy]/this.vL
end

# bulk flux
function flux!(f,u, edge, this::YSZParameters)
    uk=viewK(edge,u)
    ul=viewL(edge,u)
    f[iphi]=this.eps0*(1+this.chi)*(uk[iphi]-ul[iphi])    
    
    bp,bm=fbernoulli_pm(
        (1.0 + this.mO/this.ML*this.m_par*(1.0-this.nu))
        *(log(1-ul[iy]) - log(1-uk[iy]))
        -
        this.zA*this.e0/this.T/this.kB*(
            1.0 + this.mO/this.ML*this.m_par*(1.0-this.nu)*0.5*(uk[iy]+ul[iy])
        )*(ul[iphi] - uk[iphi])
    )
    f[iy]= (
        this.DD
        *
        (1.0 + this.mO/this.ML*this.m_par*(1.0-this.nu)*0.5*(uk[iy]+ul[iy]))
        *
        this.mO*this.m_par*(1.0-this.nu)/this.vL
        *
        (bm*uk[iy]-bp*ul[iy])
    )
end


# sources
function reaction!(f,u, node, this::YSZParameters)
    f[iphi]=-(this.e0/this.vL)*(this.zA*this.m_par*(1-this.nu)*u[iy] + this.zL) # source term for the Poisson equation, beware of the sign
    f[iy]=0
end

function bstorage!(f,u,node, this::YSZParameters)
    if  node.region==1
    if this.ARO_mode
        f[iyAs]=this.mO*this.ms_par*(1.0-this.nus)*u[iyAs]/this.areaL
        f[iyOs]=this.mO*this.OC*u[iyOs]/this.areaL
        f[iyOmins]=this.mO*this.ms_par*(1.0-this.nus)*u[iyOmins]/this.areaL
    else
        f[iyAs]=this.mO*this.ms_par*(1.0-this.nus)*u[iyAs]/this.areaL
        f[iyOs]=this.mO*this.ms_par*(1.0-this.nus)*u[iyOs]/this.areaL
        f[iyOmins]=this.mO*this.ms_par*(1.0-this.nus)*u[iyOmins]/this.areaL
    end
end

end

# surface reaction + adsorption
function breaction!(f,u,node,this::YSZParameters)
    if  node.region==1
        electroR=electroreaction(this,u)
        oxide_ads = exponential_oxide_adsorption(this, u)
        gas_ads = exponential_gas_adsorption(this, u)
        
        #
        g_ARO = similar(f)
        
        g_ARO[iy]= this.mO*oxide_ads
        g_ARO[iyAs]= - this.mO*electroR - this.mO*oxide_ads
        g_ARO[iyOs]= this.mO*electroR - this.mO*2*gas_ads
        g_ARO[iyOmins]= 0  #this.mO*(0.1 - u[iyOmins])
        g_ARO[iphi]=0

        #h_template = similar(f)
        #[@show ii,breac.value for (ii,breac) in enumerate(h_template)]
        rema = this.stoichiometric_matrix
        rates = reaction_rates(u, this)
        h_template = (this.mO*rema*rates)[1:end-1]
        #[@show ii,breac.value for (ii,breac) in enumerate(h_template)]

        if this.ARO_mode
            f = g_ARO
        else
            for (ii, hh) in enumerate( h_template)
              f[ii] = -hh
              #@show hh.value
            # f[1] = 0
            # f[2] = h_template[2]#-this.mO*u[2]
            # f[3] = h_template[3]#-this.mO*u[3]
            # f[4] = h_template[4]#-this.mO*u[4]
            # f[5] = h_template[5]#-this.mO*u[5]
          end
        end 
    else
        f[iy]=0
        f[iphi]=0
    end

end

function breaction_identity_test!(f,u,node,this::YSZParameters)
    if  node.region==1
        f[iy]=   u[iy]
        f[iyAs]= u[iyAs]
        f[iyOs]= u[iyOs]
        f[iyOmins]= u[iyOmins] 
        f[iphi]=  0
    else
        f[iy]=0
        f[iphi]=0
    end

end

# surface reactions
function exponential_oxide_adsorption(this::YSZParameters, u; debug_bool=false)
    if this.rA > 0
        # O-2(y) + V(s) => O-2(s) + V(y)
        if this.kappaA > 0
          # LoMA
          the_fac = (
                (u[iy]*(1-u[iy]))
                *
                (u[iyAs]*(1-u[iyAs]))               
              )^(this.SA*this.kappaA)
        else  
          the_fac = 1.0
        end
        rate = (
            (this.rA/this.SA)*the_fac
            *(
                exp(-this.betaA*this.SA*this.DGA/(this.kB*this.T))
                *(
                  (u[iy]/(1-u[iy]))
                  *
                  ((1-u[iyAs])/u[iyAs])
                )^(this.betaA*this.SA)
                - 
                exp((1 - this.betaA)*this.SA*this.DGA/(this.kB*this.T))
                *(
                  (u[iy]/(1-u[iy])) 
                  *
                  ((1-u[iyAs])/u[iyAs])
                )^(-(1 - this.betaA)*this.SA)
                
            )
        )
    else
      the_fac = 0
      rate=0
    end
    if debug_bool
      print("  A > ")
      a_reac = (
                  (u[iy]/(1-u[iy]))
                  *
                  ((1-u[iyAs])/u[iyAs])
                )^(this.betaA*this.SA)
      a_prod = (
                  (u[iy]/(1-u[iy]))
                  *
                  ((1-u[iyAs])/u[iyAs])
                )^(-(1 - this.betaA)*this.SA)
      @show the_fac, rate, a_reac, a_prod
    end
    return rate
end

function electroreaction(this::YSZParameters, u; debug_bool=false)
    if this.rR > 0
        # O(s) + 2e-(s) => O-2(s)
        if this.kappaR > 0
          # LoMA
          the_fac = (
               (u[iyOs]*(1-u[iyOs]))
               *
               (u[iyAs]*(1-u[iyAs]))               
            )^(this.SR*kappa.R)
        else  
          the_fac = 1
        end
        rate = (
            (this.rR/this.SR)*the_fac
            *(
                exp(-this.betaR*this.SR*this.DGR/(this.kB*this.T))
                *(u[iyAs]/(1-u[iyAs]))^(-this.betaR*this.SR)
                *(u[iyOs]/(1-u[iyOs]))^(this.betaR*this.SR)
                - 
                exp((1.0-this.betaR)*this.SR*this.DGR/(this.kB*this.T))
                *(u[iyAs]/(1-u[iyAs]))^((1.0-this.betaR)*this.SR)
                *(u[iyOs]/(1-u[iyOs]))^(-(1.0-this.betaR)*this.SR)
            )
        )
    else
      the_fac = 0
      rate = 0
    end
    if debug_bool
      print("  R > ")
      @show the_fac, rate
    end
    return rate
end

function exponential_gas_adsorption(this::YSZParameters, u; debug_bool=false)
    if this.rA > 0 && !(this.pO2 == 0)
        # O2(g) => 2O(s)
        if this.kappaO > 0
          # LoMA
          the_fac = (
               (this.pO2)
               *
               (u[iyOs]*(1-u[iyOs]))^2                                            
            )^(this.SO*this.kappaO)
        else  
          the_fac = 1.0
        end
        rate = (
            (this.rO/this.SO)*the_fac
            *(
                exp(-this.betaO*this.SO*this.DGO/(this.kB*this.T))
                *(u[iyOs]/(1-u[iyOs]))^(-2*this.betaO*this.SO)
                *(this.pO2)^(this.betaO*this.SO)
                - 
                exp((1.0-this.betaO)*this.SO*this.DGO/(this.kB*this.T))
                *(u[iyOs]/(1-u[iyOs]))^(2*(1.0-this.betaO)*this.SO)
                *(this.pO2)^(-(1.0-this.betaO)*this.SO)
            )
        )
    else
      the_fac = 0
      rate=0
    end
    if debug_bool
      print("  O > ")
      @show the_fac, rate, this.pO2
    end
    return rate
end



function reaction_template(u, 
                           this::YSZParameters,
                           r0::Float64,
                           DG::Float64,
                           beta::Float64,
                           S::Float64,
                           kappa::Float64,
                           gamma::Array{Int32,1};
                           debug="",
         )
    activities = activity(this, u)

    ybV = 1-u[iy]
    ysV = 1-sum(u[3:5])

    bv_count = gamma[2]
    sv_count = sum(gamma[3:5])
    if this.ARO_mode # || (debug != "") 
      activities[iyAs] = u[iyAs]/(1-u[iyAs])
      activities[iyOs] = u[iyOs]/(1-u[iyOs])
      if gamma[5] != 0
          ysV = 1
      #    println("Skip")
      elseif gamma[2] == 1
          ysV = 1-u[iyAs]
      #    @show gamma, "A"
      elseif gamma[6] == -1
          ysV = 1-u[iyOs]
      #    @show gamma, "O"
      else
          ysV = (1-u[iyAs])*(1-u[iyOs])
          sv_count = 1
      #    @show gamma, "R"
      end
      #println("INFO: ARO_mode shared surface vacancy override.")
    end
    if false
        for uu in u[2:5]
          @show uu
          if 0.0 < uu < 1.0
            @show u[2:5]
            error("ypsilon out of (0,1)")
          end
        end

        if 0.0 < ybV < 1.0
          @show ybV
          error("Bulk vacancies out of (0,1)")
        end

        if 0.0 < ysV < 1.0
          @show ysV
          error("Surf vacancies out of (0,1)")
        end
    end
  
    reac_activities = prod(activities.^gamma)
    #   exp(-beta*S*DG/(kB*T))
    #           *(
    #             (u[iy]/(1-u[iy]))
    #             *
    #             ((1-u[iyAs])/u[iyAs])
    #           )^(this.betaA*this.SA)

    L = exp(
              -beta*S*DG/this.kB/this.T
            )*reac_activities^(-beta*S)
    R = exp(
               (1-beta)*S*DG/this.kB/this.T
            )*reac_activities^((1-beta)*S)
    
    prefactor = prod(
                     vcat(
                         [1, 
                          u[iy], 
                          u[iyAs], 
                          u[iyOs], 
                          u[iyOmins], 
                          this.pO2
                         ].^abs.(gamma), 
                         [ybV, 
                          ysV
                         ].^abs.([bv_count, sv_count])
                     )
                )^(S*kappa)
    rate = r0/S*prefactor*(L - R)

      a_reac = (reac_activities)^(-beta*S)
      a_prod = (reac_activities)^((1 - beta)*S)
      the_fac = prefactor

    if debug != ""
        print("-------- $(debug) Reaction Template debug--------\n")
        @show gamma,bv_count,sv_count
        #@show S*kappa
        #@show activities
        #@show reac_activities
        print("  $(debug) > ") 
        @show the_fac.value, rate.value ,a_reac.value, a_prod.value
        @show rate
        println("-------o----------------")
    end
    return rate
end

function activity(this::YSZParameters, u)
    ybV = 1-u[iy]
    ysV = 1-sum(u[3:5])
    #[@show x.value for x in u]
    return [1, u[iy]/ybV, u[iyAs]/ysV, u[iyOs]/ysV, u[iyOmins]/ysV, this.pO2] 
end

function reaction_rates(u, this::YSZParameters)
    A = reaction_template(u, this,
                           this.rA,
                           this.DGA,
                           this.betaA,
                           this.SA,
                           this.kappaA,
                           this.gammaA;
                           #debug="A"
        )
    R = reaction_template(u, this,
                           this.rR,
                           this.DGR,
                           this.betaR,
                           this.SR,
                           this.kappaR,
                           this.gammaR,
                           #debug="R"
        )
    O = reaction_template(u, this,
                           this.rO,
                           this.DGO,
                           this.betaO,
                           this.SO,
                           this.kappaO,
                           this.gammaO,
                           #debug="O"
        )
    B = reaction_template(u, this,
                           this.rB,
                           this.DGB,
                           this.betaB,
                           this.SB,
                           this.kappaB,
                           this.gammaB,
                           #debug="B"
        )
    C = reaction_template(u, this,
                           this.rC,
                           this.DGC,
                           this.betaC,
                           this.SC,
                           this.kappaC,
                           this.gammaC,
                           #debug="C"
        )
    #[@show rate.value for rate in [A,R,O, B,C]]
    return [A,R,O, B, C]
end

#
# Transient part of measurement functional
#

function set_meas_and_get_tran_I_contributions(meas, U, sys, parameters, AreaEllyt, X)
    Qb= - integrate(sys,reaction!,U) # \int n^F            
    dphi_end = U[iphi, end] - U[iphi, end-1]
    dx_end = X[end] - X[end-1]
    dphiB=parameters.eps0*(1+parameters.chi)*(dphi_end/dx_end)
    Qs= (parameters.e0/parameters.areaL)*parameters.zA*U[iyAs,1]*parameters.ms_par*(1-parameters.nus) # (e0*zA*nA_s)
    meas[1]= AreaEllyt*( -Qs[1] -Qb[iphi]  -dphiB)
    return ( -AreaEllyt*Qs[1], -AreaEllyt*Qb[iphi], -AreaEllyt*dphiB)
end

#
# Steady part of measurement functional
#
function set_meas_and_get_stdy_I_contributions(meas, U, sys, parameters, AreaEllyt, X)
    meas[1] = AreaEllyt*(-2*parameters.e0*electroreaction(parameters, U[:, 1]))
    return AreaEllyt*(-2*parameters.e0*electroreaction(parameters, U[:, 1]))
end

function OCV_test(;
                print_bool=false,
                test=false, 
                verbose=false, 
                width=0.0005, dx_exp=-9,
                pO2=1.0, T=1073,
                prms_names_in=[],
                prms_values_in=[],
                 )
    model_label = "ysz_shn"
    model_symbol = eval(Symbol(model_label))
    bulk_species = model_symbol.bulk_species
    surface_species = model_symbol.surface_species
    surface_names = model_symbol.surface_names
    iphi = model_symbol.iphi
    iy = model_symbol.iy
    iyAs = model_symbol.iyAs
    iyOs = model_symbol.iyOs
    #iyOmins = model_symbol.iyOmins
    index_driving_species = model_symbol.index_driving_species

    
    # prms_in = [ rA, rR, DGA, DGR, beta, A ]

    # Geometry of the problem
    #AreaEllyt = 0.000201 * 0.6      # m^2   (geometrical area)*(1 - porosity)
    AreaEllyt = 0.018849556 * 0.7        # m^2 (geometrical area)*(1 - porosity)
    #width_Ellyt = 0.00045           # m     width of the half-cell
    #width_Ellyt = 0.0005           # m     width of the half-cell
    #
    dx_start = 10^convert(Float64,dx_exp)
    X=width*VoronoiFVM.geomspace(0.0,1.0,dx_start,1e-1)
    #println("X = ",X)
    #
    grid=VoronoiFVM.Grid(X)
    #
    
    # Physical condtions and parameters setting
    parameters=model_symbol.YSZParameters()
    

    # experimental setting
    parameters.pO2 = pO2
    parameters.T = T
    
    model_symbol.set_parameters!(parameters, prms_values_in, prms_names_in)

    # update the "computed" values in parameters
    parameters = model_symbol.YSZParameters_update!(parameters)
    if print_bool
      model_symbol.printfields(parameters)
    end
    
    physics=VoronoiFVM.Physics(
        data=parameters,
        num_species=size(bulk_species,1)+size(surface_species,1),
        storage=model_symbol.storage!,
        flux=model_symbol.flux!,
        reaction=model_symbol.reaction!,
        breaction=model_symbol.breaction!,
        bstorage=model_symbol.bstorage!
    )
    #

    #sys=VoronoiFVM.SparseSystem(grid,physics)
    sys=VoronoiFVM.DenseSystem(grid,physics)
    for idx in bulk_species
      enable_species!(sys, idx, [1])
    end
    for idx in surface_species
      enable_boundary_species!(sys, idx, [1])
    end

    # set boundary conditions
    model_symbol.set_typical_boundary_conditions!(sys, parameters)
    
    # get initial value of type unknows(sys) and initial voltage
    inival = model_symbol.get_typical_initial_conditions(sys, parameters)
    # @show inival[:,1]
    phi0 = parameters.phi_eq

  
    function breactionwrap(u)
      return parameters.mO*(parameters.stoichiometric_matrix*reaction_rates(u, parameters))[1:end-1]
    end
    xx = [parameters.phi_eq,
          parameters.y0_eq,
          parameters.yAs_eq,
          parameters.yOs_eq,
          parameters.yOmins_eq,
         ]
    brec_jac = ForwardDiff.jacobian(breactionwrap, xx)
    #[@show  brec_jac[ii,:] for ii=1:5]
    # return ForwardDiff.jacobian(breactionwrap, xx)
    #
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose
    control.tol_linear=1.0e-5
    control.tol_relative=1.0e-5
    control.tol_absolute=1.0e-4
    control.max_iterations=2000
    control.max_lureuse=0
    control.damp_initial=5.0e-3
    control.damp_growth=1.4

##### used for diplaying initial conditions vs steady state    
     #figure(111)
     plot_solution(inival, X, 10^9)
     steadystate = unknowns(sys)
     steadystate = model_symbol.get_typical_initial_conditions(sys, parameters)
     if test
         t=0.0
         tstep=1e-4
         tend=1.0
         while t<tend
             @show t
             @show breactionwrap(inival[:,1])
             #try 
                 solve!(steadystate, inival, sys, control=control, tstep=tstep)
             #catch e
             #    println(e)
                 #solve!(steadystate, inival, sys, control=control, tstep=tstep)
             #end
             # I_stdy  = set_meas_and_get_stdy_I_contributions(meas, steadystate, sys, parameters, AreaEllyt, X)
             #a = set_meas_and_get_tran_I_contributions(meas, steadystate, sys, parameters, AreaEllyt, X)
             #b = set_meas_and_get_tran_I_contributions(meas, inival, sys, parameters, AreaEllyt, X)
             #I_tran = [ (a[ii]-b[ii])/tstep for ii=1:3]
             #@show I_tran, I_stdy
             inival=steadystate
             tstep=1.3*tstep
             t+=tstep
             plot_solution(steadystate, X, 10^9)
         end
     else
         solve!(steadystate, inival, sys, control=control)
         @show breactionwrap(inival[:,1])
         @show breactionwrap(steadystate[:,1])
         plot_solution(steadystate, X, 10^9)
     end
     return
################
end 

function plot_solution(U, X, x_factor=10^9, x_label="", plotted_length= 5.0)
  point_marker_size = 5
  
  
  subplot(211)
  PyPlot.plot((x_factor)*X[:],U[iy,:],label="y")
  for (i, idx) in enumerate(surface_species)
    PyPlot.plot(0,U[idx,1],"o", markersize=point_marker_size, label=surface_names[i])
  end
  PyPlot.ylim(-0.5,1.1)
  PyPlot.xlim(-0.01*plotted_length, plotted_length)
  PyPlot.xlabel(x_label)
  PyPlot.legend(loc="best")
  PyPlot.grid()
  
  subplot(212)
  PyPlot.plot((x_factor)*X[:],U[iphi,:],label="phi (V)")
  PyPlot.xlim(-0.01*plotted_length, plotted_length)
  PyPlot.xlabel(x_label)
  PyPlot.legend(loc="best")
  PyPlot.grid()
end

end
