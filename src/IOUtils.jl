function ensure_parent_dir(path::AbstractString)
    mkpath(dirname(path))
    return path
end

function write_dataframe_csv(df::DataFrame, path::AbstractString)
    ensure_parent_dir(path)
    open(path, "w") do io
        println(io, join(string.(names(df)), ","))
        for row in eachrow(df)
            println(io, join(string.collect(row), ","))
        end
    end
    return path
end 

function save_summary_table(df::DataFrame, path::AbstractString)
    return write_dataframe_csv(df, path)
end


function save_pipeline_tables(pipeline_outputs::Dict{Symbol, Any}, path::AbstractString)
    tables = DataFrame[]
    for output in values(pipeline_outputs)
        pushs!(tables, pipeline_coeftable(output))
    end
    merged = vcat(tables...)
    return write_dataframe_csv(merged, path)
end