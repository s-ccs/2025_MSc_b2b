function make_multichannel(dat, dat_e, cfg::SimConfig; rng = MersenneTwister(cfg.seed))
    dat_e_multi = repeat(dat_e, cfg.n_channels, 1, 1)
    dat_multi = permutedims(repeat(dat, 1, cfg.n_channels), [2, 1])

    dat_e_multi . += cfg.channel_jitter_sd .* randn(rng, size(dat_e_multi)...)
    