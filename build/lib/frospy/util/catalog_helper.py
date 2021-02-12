import obspy


def select_from_cat(cat, cmts):
    cat_out = obspy.Catalog()
    for e in cat:
        ADD = False
        for ed in e.event_descriptions:
            if ed.text in cmts:
                ADD = True
        if ADD is True:
            cat_out += e
    return cat_out
