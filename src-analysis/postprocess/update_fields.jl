function update_fields!(fields, fds, clock, frame)

    # Set current state
    set!(fields.u, fds.u[frame])
    set!(fields.v, fds.v[frame])
    set!(fields.w, fds.w[frame])
    set!(fields.b, fds.b[frame])
    set!(fields.pNHS, fds.pNHS[frame])

    # Set next state
    l = length(fds.u)
    f2 = min(frame + 1, l)
    set!(fields.u_next, fds.u[f2])
    set!(fields.v_next, fds.v[f2])
    set!(fields.w_next, fds.w[f2])
    set!(fields.b_next, fds.b[f2])
    set!(fields.pNHS_next, fds.pNHS[f2])

    compute_background!(fields.U, fields.V, fields.W, clock)

    return nothing
end
