#!/usr/bin/env python3
"""
Fix forces.c fallback warnings to print only once per unique combination.
"""

with open('TREE/src/forces.c', 'r') as f:
    content = f.read()

# Fix bond parameter lookup - add static warning arrays
old_bond = '''    /* Fall back to single-bond entry if double/triple not found */
    if (order > 1) {
        for (int i = 0; i < BOND_TABLE_LEN; i++) {
            const BondParam *p = &BOND_TABLE[i];
            if (p->Za == Za && p->Zb == Zb && p->order == 1) {
                fprintf(stderr, "forces: WARNING (audit fix F3): no BOND_TABLE entry for Z%d-Z%d order %d - falling back to single-bond parameters instead of the requested bond order\\n", Za, Zb, order);
                *out = *p;
                return 1;
            }
        }
    }

    /* Geometric fallback: use sum of covalent radii, generic k */
    const Element *ea = pt_element(Za);
    const Element *eb = pt_element(Zb);
    if (ea && eb) {
        fprintf(stderr, "forces: WARNING (audit fix F3): no BOND_TABLE entry for %s(Z%d)-%s(Z%d) order %d - geometric fallback: r0 from covalent-radii sum, generic k = 20 eV/A^2\\n", ea->symbol, Za, eb->symbol, Zb, order);
        out->Za = Za; out->Zb = Zb; out->order = order;
        out->r0 = ea->covalent_radius + eb->covalent_radius;
        out->k  = 20.0;  /\\* generic, eV/A\\u00b2 *\\/
        return 1;
    }'''

new_bond = '''    /* Fall back to single-bond entry if double/triple not found */
    if (order > 1) {
        for (int i = 0; i < BOND_TABLE_LEN; i++) {
            const BondParam *p = &BOND_TABLE[i];
            if (p->Za == Za && p->Zb == Zb && p->order == 1) {
                static volatile int bond_order_warned[118][118][4] = {{0}};
                if (!bond_order_warned[Za][Zb][order]) {
                    bond_order_warned[Za][Zb][order] = 1;
                    fprintf(stderr, "forces: WARNING (audit fix F3): no BOND_TABLE entry for Z%d-Z%d order %d - falling back to single-bond parameters instead of the requested bond order\\n", Za, Zb, order);
                }
                *out = *p;
                return 1;
            }
        }
    }

    /* Geometric fallback: use sum of covalent radii, generic k */
    const Element *ea = pt_element(Za);
    const Element *eb = pt_element(Zb);
    if (ea && eb) {
        static volatile int bond_geom_warned[118][118][4] = {{0}};
        if (!bond_geom_warned[Za][Zb][order]) {
            bond_geom_warned[Za][Zb][order] = 1;
            fprintf(stderr, "forces: WARNING (audit fix F3): no BOND_TABLE entry for %s(Z%d)-%s(Z%d) order %d - geometric fallback: r0 from covalent-radii sum, generic k = 20 eV/A^2\\n", ea->symbol, Za, eb->symbol, Zb, order);
        }
        out->Za = Za; out->Zb = Zb; out->order = order;
        out->r0 = ea->covalent_radius + eb->covalent_radius;
        out->k  = 20.0;  /\\* generic, eV/A\\u00b2 *\\/
        return 1;'''

content = content.replace(old_bond, new_bond)

# Fix angle parameter lookup - add static warning array
old_angle = '''    /* Generic tetrahedral fallback - k matches the corrected table scale
     * (2x AMBER_parm x KCAL_MOL_TO_EV convention, see table comment above).
     * An earlier version left this at the OLD pre-correction value (0.60)
     * after the table itself was fixed, silently giving any untabulated
     * angle type inconsistent, too-soft physics relative to every
     * tabulated entry - using a representative generic AMBER value
     * (40 kcal/mol/rad^2, the CT-CT-CT constant) here instead. */
    fprintf(stderr, "forces: WARNING (audit fix F3): no ANGLE_TABLE entry for angle Z%d-Z%d-Z%d - generic tetrahedral fallback (109.47 deg, k = 3.469 eV/rad^2)\\n", Za, Zb, Zc);
    out->Za = Za; out->Zb = Zb; out->Zc = Zc;
    out->theta0 = DEG2RAD(109.47);
    out->k      = 2.0 * 40.0 * KCAL_MOL_TO_EV;  /\\* = 3.469 eV/rad^2 *\\/
    return 0;'''

new_angle = '''    /* Generic tetrahedral fallback - k matches the corrected table scale
     * (2x AMBER_parm x KCAL_MOL_TO_EV convention, see table comment above).
     * An earlier version left this at the OLD pre-correction value (0.60)
     * after the table itself was fixed, silently giving any untabulated
     * angle type inconsistent, too-soft physics relative to every
     * tabulated entry - using a representative generic AMBER value
     * (40 kcal/mol/rad^2, the CT-CT-CT constant) here instead. */
    static volatile int angle_warned[118][118][118] = {{0}};
    if (!angle_warned[Za][Zb][Zc]) {
        angle_warned[Za][Zb][Zc] = 1;
        fprintf(stderr, "forces: WARNING (audit fix F3): no ANGLE_TABLE entry for angle Z%d-Z%d-Z%d - generic tetrahedral fallback (109.47 deg, k = 3.469 eV/rad^2)\\n", Za, Zb, Zc);
    }
    out->Za = Za; out->Zb = Zb; out->Zc = Zc;
    out->theta0 = DEG2RAD(109.47);
    out->k      = 2.0 * 40.0 * KCAL_MOL_TO_EV;  /\\* = 3.469 eV/rad^2 *\\/
    return 0;'''

content = content.replace(old_angle, new_angle)

with open('TREE/src/forces.c', 'w') as f:
    f.write(content)

print("Forces.c fallback warnings fixed!")