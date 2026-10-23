# Integration checks on the configured ensemble; no simulation is rerun.
using Test
include(joinpath(@__DIR__,"..","src","figures","cross_so4_thermodynamics.jl"))
dirs,set_name=resolve_run_set_args(["--set","full_ensemble_2026_08_01"];min_dirs=2)
selection,extremes=process_representatives(dirs,set_name)
@testset "Outcome-only case selection" begin
  @test count(selection.selected)==2
  for regime in REGIMES[1:2]
    sub=selection[selection.regime .== regime,:]
    @test only(sub.distance_from_median[sub.selected])==minimum(sub.distance_from_median)
  end
end

rep=first(eachrow(selection[selection.selected,:]))
run=load_process_run(last(extremes),rep.seed)
ps=run.ps
@testset "Substrate-matched route accounting" begin
  mapping=matched_vfa_routes(ps)
  @test mapping[7].substrate=="Butyrate" && mapping[7].route=="SR"
  @test mapping[15].route=="acetate_forming" && mapping[16].route=="other"
  @test mapping[8].substrate=="Propionate" && mapping[17].route=="acetate_forming"
  @test mapping[18].route=="other" && mapping[13].route=="other"
  state=copy(run.C[end,:])
  state[(ps.N+1):(ps.N+ps.M)].=0.001
  C=permutedims(hcat(state,state))
  routes=integrate_vfa_routes(C,[0.,10.],ps;n_points=2)
  @test nrow(routes)==9
  cat_of=Dict(rid=>cat for (cat,rids) in CAT for rid in rids)
  on,off=same_state_category_rates(state,ps,cat_of)
  @test sum(routes.actual_substrate_mol[routes.route.=="SR"]) ≈ 10*on["VFA_SR_to_acet"]
  @test sum(routes.reference_substrate_mol[routes.route.=="SR"]) ≈ 10*off["VFA_SR_to_acet"]
  for g in groupby(routes,:substrate)
    @test sum(g.substrate_fraction) ≈ 1
  end
  @test all(x->isnan(x)||0<=x<=1,routes.inhibition)
  @test all(x->isnan(x)||isapprox(x,1;atol=1e-12),routes.effective_dG_reference_coverage)
end
@testset "ODE-consistent pathway budget and exact factorization" begin
  for ti in (2,length(run.T)÷2,length(run.T))
    state=copy(run.C[ti,:])
    rates=zeros(ps.N,ps.O)
    dx=zeros(length(state))
    Analysis.AssemblyCook.full_dynamics!(dx,state,ps,rates,run.T[ti],true)
    budget=pathway_budget_state(state,ps)
    for (ci,cat) in enumerate(CAT_ORDER)
      expected=sum(rates[i,rid]*max(state[i],0)/Analysis.NA for i in 1:ps.N for rid in CAT[cat])
      @test isapprox(budget[ci,7],expected;rtol=1e-12,atol=1e-25)
    end
  end
  state=copy(run.C[end,:])
  state[1:ps.N].=100.0
  state[(ps.N+ps.M+1):(2*ps.N+ps.M)].=0.33
  @test all(pathway_budget_state(state,ps)[:,5:7].==0)
  state[(ps.N+ps.M+1):(2*ps.N+ps.M)].=1.0
  state[(ps.N+1):(ps.N+ps.M)].=0.001
  C=permutedims(hcat(state,state))
  metrics,trace=integrated_pathway_budget(C,[0.,10.],ps;n_points=2)
  @test nrow(metrics)==length(CAT_ORDER)
  @test nrow(trace)==2*length(MECHANISM_CATS)
  for r in eachrow(metrics)
    isfinite(r.factorized_actual_flux_mol) || continue
    @test isapprox(r.factorized_actual_flux_mol,r.actual_flux_mol;rtol=1e-12,atol=1e-25)
    @test isapprox(r.reference_flux_mol-r.thermal_loss_mol,r.actual_flux_mol;rtol=1e-10,atol=1e-25)
  end
end
strain=findfirst(m->22 in m.Reacs,ps.mics)
j=findfirst(==(22),ps.mics[strain].Reacs)
mic=ps.mics[strain]
rxn=ps.reacs[22]
@testset "Effective driving force and direct kinetic suppression" begin
  for ti in (2,length(run.T)÷2,length(run.T))
    state=copy(run.C[ti,:])
    obs=process_observation(state,ps,strain,j)
    if isfinite(obs.theta_identity_error)
      @test obs.theta_identity_error<1e-10
    end
    if obs.flux_same_state_theta0_mol_s>0
      S=state[ps.N .+ rxn.Rct]
      limiting=minimum(S ./ rxn.stoichR)
      ratio=(1-obs.theta)*(mic.KS[j]+limiting)/(mic.KS[j]+limiting*(1+mic.kr[j]*obs.theta))
      @test obs.flux_on_mol_s/obs.flux_same_state_theta0_mol_s ≈ ratio
    end
  end
  state=copy(run.C[end,:])
  state[(ps.N+1):(ps.N+ps.M)].=0.001
  state[ps.N+first(rxn.Rct)]=0.0
  obs=process_observation(state,ps,strain,j)
  @test obs.theta==1.0
  @test obs.flux_on_mol_s==0.0
  @test isnan(obs.effective_delta_g_kj_mol)
  @test obs.boundary_status=="nonpositive_substrate"
  state[(ps.N+1):(ps.N+ps.M)].=0.001
  state[ps.N+14]=0.0
  obs=process_observation(state,ps,strain,j)
  @test obs.theta==0.0
  @test obs.flux_on_mol_s==obs.flux_same_state_theta0_mol_s
  @test isnan(obs.effective_delta_g_kj_mol)
  @test obs.boundary_status=="nonpositive_product"
end
