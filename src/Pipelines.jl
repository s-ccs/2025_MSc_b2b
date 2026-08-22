include("pipelines/standard_decoding.jl")
include("pipelines/rerp_decoding.jl")
include("pipelines/plain_b2b.jl")
include("pipelines/one_step_b2b.jl")
include("pipelines/two_step_b2b.jl")


function run_debug_pipelines(
    cfg,
    cases;
    target::Symbol = :continuous,
    cross_val_reps::Int = 3, 
)
    standard = run_standard_decoding(cfg, cases;target = target,)
    rerp = run_rerp_decoding(cfg, cases; target = target,)
    plain_b2b = run_plain_b2b(cfg, cases; target = target,)

    return (
        standard = standard,
        rerp = rerp,
        plain_b2b = plain_b2b,
    )
end 