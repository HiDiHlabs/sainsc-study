using CSV
using DataFrames
using Sainsc

stereoseq_file = ARGS[1]
signature_file = ARGS[2]
gene_file = ARGS[3]

counts = readstereoseq(stereoseq_file)

genes = String[]
open(gene_file, "r") do f
    for ln in eachline(f)
        push!(genes, ln)
    end
end

signatures = CSV.read(signature_file, DataFrame, transpose=true)
celltypes = signatures.Column1
select!(signatures, Not(:Column1))
signatures = signatures[:, names(signatures).∈[keys(counts)]]
signatures = signatures[:, names(signatures).∈[genes]]

kernel = gaussiankernel(8, 2)
kernel = convert.(Float32, kernel)

ctmap, cosine, assignment = assigncelltype(counts, signatures, kernel; celltypes=celltypes, log=true)

println("Done")
