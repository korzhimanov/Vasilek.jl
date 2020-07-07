module Limiters

"""
    van_leer(r)

Van Leer flux limiter [Van Leer, J. Comput. Phys., 14 (4), 361 (1974)].
"""
function van_leer(r)
    return (r + abs(r))/(1.0 + abs(r))
end

function generate_limiter(limiter_name)
    if limiter_name == :VanLeer
        return van_leer
    end
end

end  # module
