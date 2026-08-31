using Vasilek
using Vasilek.Collisions: collide!, workspace as collision_workspace

include(joinpath(@__DIR__, "scheme_cases.jl"))

# Reentrancy: the property the 0.2 refactor existed for, and the only one that
# had no test.
#
# The case for turning the closure factories into scheme values was, in order of
# weight, that a 2D2P step sweeps O(N²) independent lines per direction and that
# work wants `Threads.@threads` -- which the old design could not give, because
# `PFCNonUniform` captured one `f_tmp` shared by every line and `SemiLagrangian`
# held its buffer in the closure. Immutable scheme values plus an explicit
# per-task workspace are what fixed that. Nothing asserted it.
#
# The assertion is bit-for-bit against the serial result, not approximate: these
# are independent lines, so threading must not perturb a single digit. A shared
# mutable buffer shows up as a mismatch on some lines and not others, and as a
# *different* set of lines from run to run.
#
# Workspaces are allocated inside the task rather than indexed by
# `Threads.threadid()`. Tasks migrate between threads, and `maxthreadid()`
# exceeds `nthreads()` whenever an interactive pool exists -- 8 against 4 on the
# development machine -- so a `threadid()`-indexed pool is both a correctness bug
# and an out-of-bounds read waiting to happen. This is the pattern the 2D2P
# sweeps should copy.
const THREAD_LINES = 512
const THREAD_N = 256

"Run `advect!` over many independent lines, one workspace per task."
function threaded_advect(scheme, data, c)
    out = [similar(d) for d in data]
    chunks = collect(Iterators.partition(eachindex(data),
                                         cld(length(data), max(Threads.nthreads(), 1))))
    Threads.@threads for k in eachindex(chunks)
        ws = workspace(scheme, THREAD_N)          # per task, never shared
        for l in chunks[k]
            advect!(out[l], data[l], scheme, c, ws)
        end
    end
    return out
end

serial_advect(scheme, data, c) =
    [(o = similar(d); advect!(o, d, scheme, c, workspace(scheme, THREAD_N)); o) for d in data]

@testset "Threaded sweeps reproduce the serial result" begin
    if Threads.nthreads() == 1
        @info "threading test running on a single thread; it proves nothing here. " *
              "Run with `julia -t 4` (CI does) for it to mean anything."
    end
    println("  nthreads = ", Threads.nthreads(), ", maxthreadid = ", Threads.maxthreadid())

    data = [[1.0 + 0.5*sin(2π*(i-1)/THREAD_N + 0.01*l) for i = 1:THREAD_N]
            for l = 1:THREAD_LINES]

    cases = vcat(uniform_schemes(fmin = 0.0, fmax = 2.0),
                 [("PFCNonUniform",
                   PFCNonUniform(fill(0.05, THREAD_N); fmin = 0.0, fmax = 2.0))])

    for (name, scheme) in cases
        c = name == "PFCNonUniform" ? 0.02 : 0.4
        # Both directions: the two branches touch different code, and a race in
        # one would not show in the other.
        for cc in (c, -c)
            threaded = threaded_advect(scheme, data, cc)
            serial = serial_advect(scheme, data, cc)
            mismatched = count(!=(0), [Int(threaded[l] != serial[l]) for l = 1:THREAD_LINES])
            mismatched == 0 ||
                println("  MISMATCH ", name, " at c = ", cc, ": ", mismatched, " lines differ")
            @test threaded == serial
        end
    end
end

@testset "Threaded collisions reproduce the serial result" begin
    # `collide!` has the same shape and the same claim, so it gets the same test.
    v = collect(range(-6, 6; length = 129))
    data = [[exp(-(x - 0.3*sin(l))^2) for x in v] for l = 1:128]
    op = BGK(1e-1)

    serial = [(o = similar(d); collide!(o, d, op, v, 0.1, collision_workspace(op, length(v))); o)
              for d in data]

    out = [similar(d) for d in data]
    chunks = collect(Iterators.partition(eachindex(data),
                                         cld(length(data), max(Threads.nthreads(), 1))))
    Threads.@threads for k in eachindex(chunks)
        ws = collision_workspace(op, length(v))
        for l in chunks[k]
            collide!(out[l], data[l], op, v, 0.1, ws)
        end
    end
    @test out == serial
end
