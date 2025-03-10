def simulate(dt=0.1 * {{COURANT}} * {{HFAC}} * {{DR}} / {{CS}},
             Ma={{MA}}, m={{M}}, x0={{X0}}, p0={{P0}}, F={{F}}, tmax={{T}},
             H={{H}}, L={{L}}, refd={{REFD}}, cs={{CS}}):
    rho0 = (1 + Ma**2) * refd
    m_fluid = refd * L * H
    p = p0 + cs**2 * (rho0 - refd)
    x = [x0]
    dxdt = 0
    ddxddt = 0
    t = [0]
    f = [p * H]
    e = [0]
    while (t[-1] <= tmax):
        rho = m_fluid / (H * (L + x[-1]))
        p = p0 + cs**2 * (rho - refd)
        f.append(p * H)
        ddxddt = (f[-1] - F) / m
        u = dxdt + 0.5 * dt * ddxddt
        x.append(x[-1] + dt * u)
        e.append(0.5 * m * u**2 + F * (x[-1] - x0))
        dxdt += dt * ddxddt
        t.append(t[-1] + dt)
    return t, x, f, e
