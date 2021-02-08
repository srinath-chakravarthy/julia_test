# McMeeking fracture paper
# Dendritic cracking in solid electrolytes driven by lithium insertion (2020)
# Journal of Power Sources
# Markus Klinsmann , Felix E. Hildebrand, Markus Ganser,Robert M. McMeeking
using CSV, DataFrames, Statistics
using Dierckx
using Plots
using LaTeXStrings
import PhysicalConstants.CODATA2018: N_A, e, R
using Unitful
using UnitfulRecipes
## adf ---- Compute Voltages ---
# p = plot(xlabel= "Current",
#         ylabel = "Bond Stregnth",
#         xlim=(0u"mA/cm^2",5u"mA/cm^2"))
p = plot()
Faraday = e*N_A
temperature = 298u"K"
den = R * temperature
rct = [10u"Ω*cm^2", 100u"Ω*cm^2"]
i0 = zeros(length(rct))u"A/m^2"
i0 .= 2.0*den./(Faraday .* rct)
for i00 in i0
    current = collect(0.0u"A/m^2":0.1u"A/m^2":50u"A/m^2") # A/m^2 = 10 mA/cm^2
    cc = zeros(length(current))

    cc .= -current./(i00)
    eta = zeros(length(current))u"V"
    eta .= 2.0 * den/Faraday .* asinh.(cc)
    sigma_bond = collect(0u"N/m^2":10e6u"N/m^2":10000e6u"N/m^2")
    omega = 13e-6u"m^3/mol"
    eta_c = zeros(length(sigma_bond))u"V"
    eta_c .= -sigma_bond .* omega ./ Faraday
    sigma_crit = zeros(length(current))u"N/m^2"
    eta_crit = ones(length(current))u"V"
    eta_crit .= eta_crit .* 10.0
    for (i,e) in enumerate(eta)
        if i > 1
            for (j,ec) in enumerate(eta_c)
                if (j > 1)
                    if (e - ec) <= 0.0u"V"
                        continue
                    else
                        if (i == 2)
                            println((e-ec))
                        end
                        sigma_crit[i] = sigma_bond[j]
                        eta_crit[i] = ec
                        break;
                    end
                end
            end
        end
    end
    println(i00/2.0)
    plot!(p , current .|> u"mA/cm^2", sigma_crit .|> u"MPa")
end

png(p, "bond_strength_limits.png")
# for (i, c) in enumerate(current)
#     e = asinh(ustrip(Float64, u"A/m^2", c))u"Unit
#     eta[i] = den/Faraday * e
# end
# # eta .= 0.5*(den/Faraday) .* asinh.(-current)
