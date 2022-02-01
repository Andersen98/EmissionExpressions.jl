module EmissionExpressions


using Symbolics

export parameters, spectrum_variable, f_actual_experimental_matlab_fitting, f_experimental_paper_AQT
export f_CUI_background_paper

parameters = @variables γ_QD γ_SP ω_QD ω_SP Ω_R
spectrum_variable = Symbolics.variable(:ω)


ex_actual_experimental_matlab_fitting() = let (a,b,c,d,G) = parameters, x_ev = spectrum_variable
    A = im*(x_ev-c + (c-d));
    B = (a+b)/4 - im*(c-d)/2 - im*(x_ev-c);
#    G = sqrt((e/2)^2 + ((c-d)/2)^2 - ((a-b)/4)^2);
    a/π * abs2((b/2-A)/(B^2+G^2))
end

f_expr_actual_experimental_matlab_fitting() = build_function(ex_actual_experimental_matlab_fitting(),γ_QD,γ_SP,ω_QD,ω_SP,Ω_R,spectrum_variable)


"""
    f_actual_experimental_matlab_fitting() -> f(γ_QD,γ_SP,ω_QD,ω_SP,Ω_R,ω)

Returns a function representing the fit code in https://github.com/mollyamay/strong-coupling-modeling/blob/b95b4dde3cbabd8a8d1dfe4266dbfb9a2ab0683c/f.m

In this expression, we assume that `G` corresponds to Ω_R (see source code)
Also, we removed a scaling parameter and an offset parameter used in the original fit code. 
This means that we are free to offset and scale as necessary.
"""
f_actual_experimental_matlab_fitting() = eval(f_expr_actual_experimental_matlab_fitting())




ex_experimental_paper_AQT() = let (γ_QD,γ_SP,ω_QD,ω_SP,Ω_R) = parameters, ω = spectrum_variable
    δ = ω_SP - ω_QD;
    γ_QD/(2π) * abs2( (γ_SP/2 - im*(ω - ω_SP)) / (((γ_SP + γ_QD)/4 + im*(δ)/2 - im*(ω-ω_QD))^2 + Ω_R^2) )
end

f_expr_experimental_paper_AQT() = build_function(ex_experimental_paper_AQT(),γ_QD,γ_SP,ω_QD,ω_SP,Ω_R,spectrum_variable)

"""
    f_experimental_paper_AQT() -> f(γ_QD,γ_SP,ω_QD,ω_SP,Ω_R,ω)
Returns function representing equation 12 of:
May, M. A., Fialkow, D., Wu, T., Park, K., Leng, H., Kropp, J. A., Gougousi, T., Lalanne, P., Pelton, M., & Raschke, M. B. (2020). Nano‐Cavity QED with Tunable Nano‐Tip Interaction. Advanced Quantum Technologies, 3(2), 1900087. https://doi.org/10.1002/qute.201900087
"""
f_experimental_paper_AQT() = eval(f_expr_experimental_paper_AQT())

#in the paper γ = γ_QD/2, κ =γ_SP/2 are defined as half of the lifetimes
ex_CUI_background_paper() = let γ = γ_QD/2, κ = γ_SP/2, ω_0 = ω_QD, ω_c = ω_SP, g = Ω_R, ω = spectrum_variable
    Δ = ω_0 -ω_c;
    K = κ + γ;
    Γ = κ - γ;
    γ/π * abs2( (κ - im*(ω - ω_c))/ (  (K/2 - im*(2ω - ω_0 - ω_c)/2)^2 +g^2    )   )
end


f_expr_CUI_background_paper() = build_function(ex_CUI_background_paper(),γ_QD,γ_SP,ω_QD,ω_SP,Ω_R,spectrum_variable)

"""
    f_CUI_background_paper() -> f(γ_QD,γ_SP,ω_QD,ω_SP,Ω_R,ω)

Returns function representing equation 14 of:
Cui, G., & Raymer, M. G. (2006). Emission spectra and quantum efficiency of single-photon sources in the cavity-QED strong-coupling regime. Physical Review A, 73(5), 053807. https://doi.org/10.1103/PhysRevA.73.053807

"""
f_CUI_background_paper() = eval(f_expr_CUI_background_paper())


end

