using PlutoUI
using Markdown

function simulation_controls()
    PlutoUI.combine() do Child

        md"""
        ### Simulation parameters
        
        **Number of trials:**
        $(Child(
            :n_trials,
            PlutoUI.Slider(
                100:100:2000;
                default = 400,
                show_value = true,
            ), 
        ))


        **Mean inter-event interval:**
        $(Child(
            :overlap_interval_ms,
            PlutoUI.Slider(
                150.0:25.0:800.0;
                show_value = true,
            ),
        ))

        **Onset condition bias:**
        $(Child(
            :onset_condition_bias,
            PlutoUI.Slider(
                -1.0:0.05:1.0;
                default = -0.6,
                show_value = true,
            ),
        ))

        **Confound strength (ρ):**
        $(Child(
            :rho,
            PlutoUI.Slider(
                0.0:0.05:1.0;
                default = 0.8,
                show_value = true,
            ),
        ))

        **EEG noise level (σ):**
        $(Child(
            :noiselevel,
            PlutoUI.Slider(
                0.0:0.05:1.0;
                default = 0.2,
                show_value = true,
            ),
        ))

        **Channel noise:**
        $(Child(
            :channel_noise_sd,
            PlutoUI.Slider(
                0.0:0.02:1.0;
                default = 0.2,
                show_value = true,
            ),
        ))

        
        **N170 intercept (β₀_n170):**
        $(Child(
            :β0_n170,
            PlutoUI.Slider(
                0.0:0.05:10.0;
                default = 5.0,
                show_value = true,
            ),
        ))

        **P300 intercept (β₀_p300):**
        $(Child(
            :β0_p300,
            PlutoUI.Slider(
                0.0:0.05:10.0;
                default = 5.0,
                show_value = true,
            ),
        ))

        **Condition effect (β_condition):**
        $(Child(
            :β_condition,
            PlutoUI.Slider(
                0.0:0.05:10.0;
                default = 3.0,
                show_value = true,
            ),
        ))

        **Continuous effect (β_continuous):**
        $(Child(
            :β_continuous,
            PlutoUI.Slider(
                0.0:0.1:6.0;
                default = 1.0,
                show_value = true,
            ),
        ))

        **Shift onset by one:**
        $(Child(
            :shift_onset,
            PlutoUI.CheckBox(default =  true),
        ))
        """
    end
end