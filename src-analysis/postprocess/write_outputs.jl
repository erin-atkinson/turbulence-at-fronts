function write_outputs(output_fts, output_fields, iteration, t)
    for (fts, field) in zip(output_fts, output_fields)
        set!(fts, field, iteration, t)
    end
    return nothing
end
