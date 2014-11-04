module Burgers

const xmin = 0.0
const xmax = 4.0

function u_initial(nx)
    dx = (xmax - xmin) / (nx - 1)
    u = Array(Float64, nx)
    for i in 1:nx
        x = xmin + (i-1)*dx
        u[i] = x<2 ? 1 : 0
    end
    return u
end

computeF(u) = u .^ 2 / 2

function maccormack(u, nt, dt, dx, epsilon)
    nx = length(u)
    un = Array(Float64, nx,nt)
    un[:,1] = u[:]
    for n in 2:nt
        # Predictor
        F = computeF(u)
        ustar = similar(u)
        ustar[1] = u[1]
        for i in 2:nx-1
            ustar[i] = (u[i] - dt/dx * (F[i+1] - F[i])
                        + epsilon * (u[i+1] - 2*u[i] + u[i-1]))
        end
        ustar[end] = u[end]
        # Corrector
        Fstar = computeF(ustar)
        un[1,n] = u[1]
        for i in 2:nx
            un[i,n] = 1/2 * (u[i] + ustar[i]) - dt/dx * (F[i] - F[i-1])
        end
        u[:] = un[:,n]
    end
    return un
end



function main()
    nx = 81
    nt = 70
    dx = (xmax - xmin) / (nx - 1)
    u = u_initial(nx)
    sigma = 0.5
    epsilon = 0.5
    dt = sigma*dx
    un = maccormack(u,nt,dt,dx,epsilon)
    # Output solution to file
    open("burgers.dat", "w") do f
        for n=1:nt
            for i=1:nx
                println(f, "$n $i $(un[i,n])")
            end
            println(f)
        end
    end
end

main()

end
