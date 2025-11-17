function update_fields!(fields, fds, clock, frame; skip_update)

    # Set current state
    :u in skip_update || set!(fields.u, fds.u[frame])
    :v in skip_update || set!(fields.v, fds.v[frame])
    :w in skip_update || set!(fields.w, fds.w[frame])
    :b in skip_update || set!(fields.b, fds.b[frame])
    :pNHS in skip_update || set!(fields.pNHS, fds.pNHS[frame])

    # Set next state
    l = length(fds.u)
    f2 = min(frame + 1, l)
    :u_next in skip_update || set!(fields.u_next, fds.u[f2])
    :v_next in skip_update || set!(fields.v_next, fds.v[f2])
    :w_next in skip_update || set!(fields.w_next, fds.w[f2])
    :b_next in skip_update || set!(fields.b_next, fds.b[f2])
    :pNHS_next in skip_update || set!(fields.pNHS_next, fds.pNHS[f2])

    compute_background!(fields.U, fields.V, fields.W, clock)

    return nothing
end
