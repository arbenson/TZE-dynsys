function _forward_euler(f, x_curr::Vector{Float64}, h::Float64)
    return x_curr + h * f (x_curr)
end

function _RK4(f, x_curr::Vector{Float64}, h::Float64)
    k1 = h * f(x_curr)
    k2 = h * f(x_curr + k1 / 2)
    k3 = h * f(x_curr + k2 / 2)
    k4 = h * f(x_curr + k3)
    return x_curr + (k1 + k2 + k3 + k4) / 6
end

"""
forward_euler
-------------
explicit Forward Euler method

Input parameters:
- h::Float64: step size

Returns:
A function that takes as input a derivative function and a vector and returns
the next iterate.  
"""
function forward_euler(h::Foat64)
    ret(f, x_curr::Float64) = _forward_euler(f, x, h)
    return ret
end


"""
RK4
-------------
Fourth-order explicit Runge-Kutta method

Input parameters:
- h::Float64: step size

Returns:
A function that takes as input a derivative function and a vector and returns
the next iterate.  
"""
function RK4(h::Foat64)
    ret(f, x_curr::Float64) = _RK4(f, x, h)
    return ret
end
