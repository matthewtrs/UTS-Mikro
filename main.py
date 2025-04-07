# --- Start of the corrected code ---
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import math
import traceback
import numpy as np
import re # Kept import just in case

# --- Matplotlib Integration ---
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
# ---------------------------

# --- Design Constants & Parameters ---
VTHN = 0.4
VTHP = -0.4
KP_PMOS = 50e-6
KN_NMOS = 100e-6
BETA1_PLACEHOLDER = KN_NMOS

# --- Helper Function for Lambda ---
# --- Using the CORRECTED concatenation logic ---
def get_lambda(base_lambda, z_digit):
    if str(base_lambda) == "0.05": prefix = "0.05"
    elif str(base_lambda) == "0.01": prefix = "0.01"
    else: prefix = str(base_lambda) # Fallback
    lambda_str = prefix + str(z_digit)
    try:
        val = float(lambda_str)
        return max(1e-9, val) # Ensure positive
    except ValueError:
        print(f"Error: Could not convert lambda string '{lambda_str}' to float.")
        return 1e-9
# --- END CORRECTED get_lambda ---

# --- Format Helper ---
# --- UNCHANGED ---
def format_value(value, unit="", precision=4, use_sci_limits=(1e-4, 1e5)):
    sci_low, sci_high = use_sci_limits
    if isinstance(value, (int, float)):
        if not math.isfinite(value): return f"{value} {unit}"
        abs_val = abs(value)
        if abs_val == 0: fmt = f".{precision}f"
        elif abs_val < sci_low or abs_val >= sci_high: fmt = ".3e"
        else: fmt = f".{precision}f"
        return f"{value:{fmt}} {unit}"
    return f"{value} {unit}"
# --- END UNCHANGED ---

# --- LaTeX Formula Dictionary --- ## MODIFIED SECTION ##
# --- Start Replace latex_formulas dictionary ---
latex_formulas = {
    1: r"$C_c > \frac{2.2}{10} C_L$",
    2: r"$I_5 = SR \cdot C_c$",
    3: r"$\frac{W_3}{L_3}=\frac{W_4}{L_4}=\frac{I_5}{K'_p[V_{DD}-V_{in(max)}-|V_{T3}|+V_{T1}]^2}$", # Matches image
    4: r"$g_{m1} = 2\pi \cdot GB \cdot C_c$",
    5: r"$S_1 = \frac{g_{m1}^2}{K'_n I_5}$",
    # --- UPDATED Step 6 Formula ---
    6: r"$V_{DS5} = V_{in(min)} - V_{SS} - \sqrt{\frac{I_5}{\beta_1}} - V_{T1(max)}$",
    # --- UPDATED Step 7 Formula ---
    7: r"$S_5 = \frac{2 I_5}{K'_n (V_{DS5})^2}$", # Uses VDS5 from Step 6
    8: r"$g_{m6} = 2.2 g_{m2} \left(\frac{C_L}{C_c}\right)$",
    9: r"$g_{m4} = g_{m3} = \sqrt{2 K'_p S_3 I_4} \quad (I_4=I_5/2)$",
    10: r"$S_6 = S_4 \cdot \frac{g_{m6}}{g_{m4}}$",
    11: r"$I_6 = Ratio \cdot I_5 \quad (\text{Ratio=Chosen})$",
    12: r"$S_7 = S_5 \left(\frac{I_6}{I_5}\right)$", # Uses S5 from Step 7
    13: r"$A_v = \frac{g_{m2} g_{m6}}{(g_{ds2}+g_{ds4})(g_{ds6}+g_{ds7})}$"
}
# --- End Replace latex_formulas dictionary ---

# --- Calculation Function (design_opamp) --- ## MODIFIED Step 3 ##
# (Added explicit gm4 calculation step, renumbered steps)
# (Based on the version BEFORE the gm4 step was added, to ensure no logic was lost)
# --- Start Replace design_opamp function ---
# --- Calculation Function (design_opamp) --- ## MODIFIED Step 3, 6, 7 ##
# (Uses image formula for S3 and VDS5 -> S5)
# --- Start Replace design_opamp ---
def design_opamp(xyz_str, target_cc, target_sr, specs_from_gui): # Receive specs from GUI
    results = {'steps': [], 'specs': {}, 'performance': {}, 'warnings': []}
    steps_data = results['steps']
    # Use specs passed from GUI, which include ratio and method
    specs = specs_from_gui.copy()
    perf = results['performance']
    warnings = results['warnings']

    # Retrieve necessary parameters from specs dict
    i6_method = specs.get('i6_method', 'ratio') # Default to ratio if missing
    i6_i5_ratio = specs.get('i6_i5_ratio', 5.0) # Default to 5 if missing

    # Basic XYZ setup (as before)
    if len(xyz_str) != 3 or not xyz_str.isdigit(): raise ValueError("Input must be exactly 3 digits (XYZ).")
    x = int(xyz_str[0]); y = int(xyz_str[1]); z = int(xyz_str[2])

    # --- Specs (use values already in specs dict) ---
    spec_gb_hz = specs['GB']; spec_vdd = specs['VDD']; spec_vss = specs['VSS']
    spec_sr_min_req = specs['SR_min']; spec_cl = specs['CL']
    spec_lambda_p = specs['lambda_p']; spec_lambda_n = specs['lambda_n']
    spec_icmr_min = specs['icmr_min_spec']; spec_icmr_max = specs['icmr_max_spec']
    spec_av_db = specs['spec_av_db']

    # --- Placeholders & Estimates ---
    est_vin_max = spec_icmr_max; specs['Vin(max)'] = est_vin_max
    est_vin_min = spec_icmr_min; specs['Vin(min)'] = est_vin_min # Use ICMR min directly
    est_vt03_abs = abs(VTHP); specs['|VT3|'] = est_vt03_abs
    est_vt1 = VTHN; specs['VT1'] = est_vt1
    est_vt1_max = VTHN + 0.0; specs['VT1(max)'] = est_vt1_max
    est_vt03_max = abs(VTHP) + 0.0; specs['|VT03|(max)'] = est_vt03_max
    est_vt1_min = VTHN - 0.0; specs['VT1(min)'] = est_vt1_min

    calculated_values = {}

    # --- Calculation Steps with Variable Tracking ---
    try:
        # Step 1: Cc
        calc_cc = target_cc; cc_min_req = 0.22 * spec_cl
        if calc_cc < cc_min_req: warnings.append(f"Chosen Cc below min.")
        calculated_values['Cc'] = calc_cc; step_vars = {'CL': spec_cl}
        steps_data.append({'step': 1, 'Langkah': 'Determine Cc (Slider/Input)', 'Hasil': calc_cc, 'Unit': 'F', 'vars': step_vars})

        # Step 2: I5
        calc_sr = target_sr
        if calc_sr < spec_sr_min_req: warnings.append(f"Chosen SR below min.")
        calc_i5 = calc_sr * calc_cc
        calculated_values['I5'] = calc_i5
        i1=i2=i3=i4=calc_i5/2.0; calculated_values['I1']=i1; calculated_values['I3']=i3; calculated_values['I4']=i4
        step_vars = {'SR': calc_sr, 'Cc': calc_cc}
        steps_data.append({'step': 2, 'Langkah': 'Determine I5 (SR*Cc)', 'Hasil': calc_i5, 'Unit': 'A', 'vars': step_vars})

        # Step 3: S3 (Using Image Formula)
        vov3_max_term = spec_vdd - est_vin_max - est_vt03_abs + est_vt1; calc_s3 = float('inf')
        if vov3_max_term <= 0: warnings.append("Step 3: Term [VDD-Vin(max)-|VT3|+VT1] <= 0.")
        elif KP_PMOS <= 0: warnings.append("Step 3: K'p <= 0.")
        else:
            denominator_term = KP_PMOS * (vov3_max_term ** 2)
            if denominator_term == 0: warnings.append("Step 3: Denominator is zero.")
            else: calc_s3 = calc_i5 / denominator_term
        calculated_values['S3'] = calc_s3; calculated_values['S4'] = calc_s3
        calculated_values['Vov3_max_calc'] = vov3_max_term
        step_vars_s3 = {'I5': calc_i5, 'K\'p': KP_PMOS, 'VDD': spec_vdd, 'Vin(max)': est_vin_max, '|VT3|': est_vt03_abs, 'VT1': est_vt1}
        steps_data.append({'step': 3, 'Langkah': 'Determine S3=(W/L)3 (P) (Image Formula)', 'Hasil': calc_s3, 'Unit': '', 'vars': step_vars_s3})

        # Step 4: gm1
        calc_gm1 = spec_gb_hz * 2 * math.pi * calc_cc
        calculated_values['gm1'] = calc_gm1; calculated_values['gm2'] = calc_gm1
        step_vars_gm1 = {'GB': spec_gb_hz, 'Cc': calc_cc}
        steps_data.append({'step': 4, 'Langkah': 'Determine gm1=gm2', 'Hasil': calc_gm1, 'Unit': 'S', 'vars': step_vars_gm1})

        # Step 5: S1
        calc_s1 = calc_gm1**2 / (KN_NMOS * calc_i5) if KN_NMOS != 0 and calc_i5 != 0 else float('inf')
        calculated_values['S1'] = calc_s1; calculated_values['S2'] = calc_s1
        calc_vov1 = math.sqrt(calc_i5 / (KN_NMOS * calc_s1)) if KN_NMOS > 0 and calc_s1 > 0 and calc_i5 >=0 else 0
        calculated_values['Vov1'] = calc_vov1
        step_vars_s1 = {'gm1': calc_gm1, 'K\'n': KN_NMOS, 'I5': calc_i5}
        steps_data.append({'step': 5, 'Langkah': 'Determine S1=(W/L)1 (N)', 'Hasil': calc_s1, 'Unit': '', 'vars': step_vars_s1})

        # Step 6: VDS5 (Using Image Formula)
        beta1 = KN_NMOS * calc_s1; calculated_values['beta1'] = beta1; term_beta = 0.0
        if beta1 > 1e-15 and calc_i5 >= 0:
             try: term_beta = math.sqrt(calc_i5 / beta1)
             except ValueError: term_beta = float('NaN'); warnings.append("Step 6: Cannot calculate sqrt(I5/beta1).")
        elif beta1 <= 1e-15: term_beta = float('inf'); warnings.append("Step 6: beta1 is near zero.")
        calc_vds5 = est_vin_min - spec_vss - term_beta - est_vt1_max
        calculated_values['Vds5'] = calc_vds5
        step_vars_vds5 = {'Vin(min)': est_vin_min, 'VSS': spec_vss, 'I5': calc_i5, 'beta1': beta1, 'VT1(max)': est_vt1_max}
        steps_data.append({'step': 6, 'Langkah': 'Calc. VDS5 (Image Formula)', 'Hasil': calc_vds5, 'Unit': 'V', 'vars': step_vars_vds5})
        chosen_vov5 = 0.2; calculated_values['Vov5_chosen_ref'] = chosen_vov5 # Keep reference Vov5
        vgs1_max = calc_vov1 + est_vt1_max; calc_vds5_sat_req = est_vin_min - spec_vss - vgs1_max; calculated_values['Vds5_sat_req'] = calc_vds5_sat_req
        if calc_vds5 < calc_vds5_sat_req: warnings.append(f"Warning: Calculated Vds5 ({format_value(calc_vds5,'V')}) < required ({format_value(calc_vds5_sat_req,'V')}) for M1/M2 sat.")

        # Step 7: S5 (Using VDS5 from Step 6)
            # --- Start Replace Step 7 (More Focused Checks) ---
    # --- Step 7: S5 (Using VDS5 calculated in Step 6) --- ## FOCUSED CHECKS ##
        vds5_for_s5_calc = calculated_values.get('Vds5', float('NaN'))
        i5_for_s5_calc = calculated_values.get('I5', float('NaN'))

        calc_s5 = float('inf') # Default to infinity

        # Check if essential inputs are finite numbers
        inputs_finite = math.isfinite(vds5_for_s5_calc) and math.isfinite(i5_for_s5_calc)
        kn_valid = KN_NMOS > 0 # K'n must be positive

        if inputs_finite and kn_valid:
            vds5_squared = vds5_for_s5_calc**2
            # Check specifically if VDS5 was exactly zero (unlikely but possible)
            if abs(vds5_squared) < 1e-18: # Effectively VDS5 == 0
                warnings.append("Step 7: VDS5 from Step 6 is zero. Cannot calculate S5.")
            else:
                # Calculate denominator
                denominator_s5 = KN_NMOS * vds5_squared
                # Final check on denominator (redundant if Kn>0 and Vds5 != 0, but safe)
                if abs(denominator_s5) < 1e-18:
                    warnings.append("Step 7: Denominator K'n*(VDS5)^2 is effectively zero.")
                else:
                    # Perform calculation
                    calc_s5 = (2 * i5_for_s5_calc) / denominator_s5
                    # Ensure calculated size is physically meaningful (non-negative)
                    if calc_s5 < 0:
                        warnings.append(f"Step 7: Calculated S5 is negative ({calc_s5}). Check I5 sign? Setting S5 to inf.")
                        calc_s5 = float('inf') # Physical size cannot be negative
        else:
            # Add warnings if inputs were invalid
            if not inputs_finite: warnings.append(f"Step 7: Invalid input(s) VDS5={vds5_for_s5_calc} or I5={i5_for_s5_calc}")
            if not kn_valid: warnings.append(f"Step 7: Invalid K'n ({KN_NMOS})")

        calculated_values['S5'] = calc_s5 # Store the result (inf if checks failed)

        # Store reference S5 calculation using chosen Vov5 (unchanged)
        chosen_vov5 = 0.2
        calc_s5_standard = (2 * i5_for_s5_calc) / (KN_NMOS * chosen_vov5**2) if KN_NMOS != 0 and chosen_vov5 > 0 and i5_for_s5_calc >= 0 else float('inf')
        calculated_values['S5_standard_ref'] = calc_s5_standard

        # Display vars and result from image formula
        step_vars_s5 = {'I5': i5_for_s5_calc, 'K\'n': KN_NMOS, 'VDS5': vds5_for_s5_calc}
        steps_data.append({
            'step': 7,
            'Langkah': 'Determine S5=(W/L)5 (N) (Image Formula)',
            'Hasil': calc_s5, # Display the calculated S5
            'Unit': '',
            'vars': step_vars_s5
        })
        # --- End Replace Step 7 (More Focused Checks) ---

        # Step 8: gm6
        gm2_val = calculated_values['gm2']; cl_val = spec_cl; cc_val = calculated_values['Cc']
        calc_gm6 = 2.2 * gm2_val * (cl_val / cc_val) if cc_val != 0 else float('inf')
        calculated_values['gm6'] = calc_gm6
        step_vars_gm6 = {'gm2': gm2_val, 'CL': cl_val, 'Cc': cc_val}
        steps_data.append({'step': 8, 'Langkah': 'Determine gm6', 'Hasil': calc_gm6, 'Unit': 'S', 'vars': step_vars_gm6})

        # Step 9: gm4
        s3_val = calculated_values['S3']; i4_val = calculated_values['I4']; calc_gm4 = 0.0
        if KP_PMOS > 0 and s3_val > 0 and i4_val >= 0 and math.isfinite(s3_val):
             try: calc_gm4 = math.sqrt(2 * KP_PMOS * s3_val * i4_val)
             except ValueError: calc_gm4 = float('NaN')
        elif i4_val == 0: calc_gm4 = 0.0
        else: calc_gm4 = float('inf')
        calculated_values['gm4'] = calc_gm4
        step_vars_gm4 = {'K\'p': KP_PMOS, 'S3': s3_val, 'I4': i4_val}
        steps_data.append({'step': 9, 'Langkah': 'Determine gm4 = gm3', 'Hasil': calc_gm4, 'Unit': 'S', 'vars': step_vars_gm4})

        # Step 10: S6_formula (for reference/use in formula I6 method)
        s4_val = calculated_values['S4']; gm6_val = calculated_values['gm6']; gm4_val = calculated_values['gm4']
        calc_s6_formula = float('inf')
        if gm4_val != 0 and math.isfinite(gm4_val): calc_s6_formula = s4_val * (gm6_val / gm4_val)
        elif not math.isfinite(gm4_val): warnings.append("Step 10: Cannot calculate S6_formula (gm4 is NaN/inf).")
        else: warnings.append("Step 10: Cannot calculate S6_formula (gm4 is zero).")
        calculated_values['S6_formula'] = calc_s6_formula # Store formula S6
        step_vars_s6 = {'S4': s4_val, 'gm6': gm6_val, 'gm4': gm4_val}
        steps_data.append({'step': 10, 'Langkah': 'Determine S6=(W/L)6 (P) (Ratio Formula)', 'Hasil': calc_s6_formula, 'Unit': '', 'vars': step_vars_s6})

        # --- Step 11: Define or Calculate I6 --- ## MODIFIED LOGIC ##
        calc_i6 = float('NaN')
        calc_s6_final = float('NaN') # S6 value consistent with chosen I6 method
        step_11_langkah = ""
        step_11_vars = {}

        if i6_method == "ratio":
            # Calculate I6 using Ratio * I5
            calc_i6 = i6_i5_ratio * calc_i5
            calculated_values['I6'] = calc_i6
            step_11_langkah = f'Define I6 = {i6_i5_ratio:.2f}*I5 (Ratio Method)'
            step_11_vars = {'I5': calc_i5, 'Ratio': i6_i5_ratio}
            # Calculate S6 based on this I6 and gm6
            calc_s6_final = float('inf')
            if KP_PMOS > 0 and calc_i6 > 0 and math.isfinite(gm6_val):
                calc_s6_final = gm6_val**2 / (2 * KP_PMOS * calc_i6)
            calculated_values['S6'] = calc_s6_final # Store S6 consistent with ratio I6

        elif i6_method == "formula":
            # Calculate I6 using gm6^2 / (2 K'p S6_formula)
            s6_formula_val = calculated_values.get('S6_formula', float('NaN'))
            calc_i6 = float('inf')
            step_11_langkah = "Calc. I6 = gm6^2 / (2 K'p S6_formula) (Formula Method)"
            step_11_vars = {'gm6': gm6_val, 'K\'p': KP_PMOS, 'S6_formula': s6_formula_val}

            if not math.isfinite(gm6_val): warnings.append("Step 11: Cannot calculate I6 (gm6 invalid).")
            elif not math.isfinite(s6_formula_val): warnings.append("Step 11: Cannot calculate I6 (S6_formula invalid).")
            elif KP_PMOS <= 0: warnings.append("Step 11: Cannot calculate I6 (K'p <= 0).")
            elif s6_formula_val <= 0: warnings.append("Step 11: Cannot calculate I6 (S6_formula <= 0).")
            else:
                denominator_i6 = 2 * KP_PMOS * s6_formula_val
                if abs(denominator_i6) < 1e-18: warnings.append("Step 11: Denominator for I6 calc is zero.")
                else: calc_i6 = gm6_val**2 / denominator_i6

            calculated_values['I6'] = calc_i6
            # In this method, the S6 used IS S6_formula
            calculated_values['S6'] = s6_formula_val
            calc_s6_final = s6_formula_val # For Vov6 calc consistency

        else: # Should not happen
            warnings.append("Step 11: Unknown I6 calculation method.")
            calculated_values['I6'] = float('NaN')
            calculated_values['S6'] = float('NaN')

        steps_data.append({'step': 11, 'Langkah': step_11_langkah, 'Hasil': calc_i6, 'Unit': 'A', 'vars': step_11_vars})
        # --- END MODIFIED STEP 11 ---

        # Step 12 (Old 11): S7
        s5_to_use = calculated_values['S5']; i6_val = calculated_values['I6']; i5_val = calculated_values['I5']
        calc_s7 = s5_to_use * (i6_val / i5_val) if i5_val != 0 and math.isfinite(i6_val) and math.isfinite(s5_to_use) else float('inf')
        calculated_values['S7'] = calc_s7
        step_vars_s7 = {'S5': s5_to_use, 'I6': i6_val, 'I5': i5_val}
        steps_data.append({'step': 12, 'Langkah': 'Determine S7=(W/L)7 (N)', 'Hasil': calc_s7, 'Unit': '', 'vars': step_vars_s7})

        # Step 13 (Old 12): Gain
        lambda_n = specs['lambda_n']; lambda_p = specs['lambda_p']
        gds2 = i1 * lambda_n; gds4 = i3 * lambda_p; gds6 = calc_i6 * lambda_p; gds7 = calc_i6 * lambda_n
        # Add checks for infinite/NaN gds before calculating ro
        gds1_sum = gds2 + gds4 if math.isfinite(gds2) and math.isfinite(gds4) else float('inf')
        gds2_sum = gds6 + gds7 if math.isfinite(gds6) and math.isfinite(gds7) else float('inf')
        ro1 = 1.0 / gds1_sum if gds1_sum > 1e-18 else float('inf')
        ro2 = 1.0 / gds2_sum if gds2_sum > 1e-18 else float('inf')
        gm2_val = calculated_values['gm2']; gm6_val = calculated_values['gm6']
        # Check for infinite/NaN before multiplying
        av1 = gm2_val * ro1 if math.isfinite(gm2_val) and math.isfinite(ro1) else float('inf')
        av2 = gm6_val * ro2 if math.isfinite(gm6_val) and math.isfinite(ro2) else float('inf')
        calc_av_lin = av1 * av2 if math.isfinite(av1) and math.isfinite(av2) else float('inf')
        calculated_values['Av_lin'] = calc_av_lin
        step_vars_gain = {'gm2': gm2_val, 'gm6': gm6_val, 'I5': i5_val, 'lambda2': lambda_n, 'lambda4': lambda_p, 'I6': i6_val, 'lambda6': lambda_p, 'lambda7': lambda_n}
        steps_data.append({'step': 13, 'Langkah': 'Determine Total Gain Av (Lin)', 'Hasil': calc_av_lin, 'Unit': '', 'vars': step_vars_gain})

    except Exception as e:
        warnings.append(f"Calculation Error: {e}"); print(f"Error: {e}"); traceback.print_exc()

    # --- Performance Summary & Warnings --- ## MODIFIED SECTION ##
    # Use final I6 and S6 values consistent with chosen method
    calc_av_db = 20 * np.log10(abs(calc_av_lin)) if calc_av_lin != 0 and math.isfinite(calc_av_lin) else -float('inf')
    perf['calc_av_lin'] = calc_av_lin; perf['calc_av_db'] = calc_av_db
    calc_power = (calculated_values.get('I5', 0) + calculated_values.get('I6', 0)) * spec_vdd if math.isfinite(calculated_values.get('I5', 0)) and math.isfinite(calculated_values.get('I6', 0)) else float('inf')
    perf['calc_power'] = calc_power

    # Vov3 Actual Calculation
    vov3_actual = float('NaN')
    s3_val = calculated_values.get('S3', float('NaN')); i3_val = calculated_values.get('I3', float('NaN'))
    if KP_PMOS > 0 and math.isfinite(s3_val) and s3_val > 0 and math.isfinite(i3_val) and i3_val >= 0:
        try: vov3_actual = math.sqrt(2 * i3_val / (KP_PMOS * s3_val))
        except ValueError: pass
    perf['calc_vov3_actual'] = vov3_actual

    # Vov5 Actual Calculation
    vov5_actual = float('NaN')
    s5_val = calculated_values.get('S5', float('NaN')); i5_val = calculated_values.get('I5', float('NaN'))
    if KN_NMOS > 0 and math.isfinite(s5_val) and s5_val > 0 and math.isfinite(i5_val) and i5_val >= 0:
         try: vov5_actual = math.sqrt(2 * i5_val / (KN_NMOS * s5_val))
         except ValueError: pass
    perf['calc_vov5_actual'] = vov5_actual

    # ICMR Calculations
    icmr_max_calc = spec_vdd - abs(vov3_actual) - abs(est_vt03_max) if math.isfinite(vov3_actual) else -float('inf')
    icmr_min_calc = spec_vss + vov5_actual + (est_vt1_min + calculated_values.get('Vov1', 0)) if math.isfinite(vov5_actual) else float('inf')
    perf['icmr_max_calc'] = icmr_max_calc; perf['icmr_min_calc'] = icmr_min_calc
    perf['calc_vov1'] = calculated_values.get('Vov1', float('NaN'))

    # Vov6 calculation using the S6 value consistent with the chosen I6 method
    s6_final_val = calculated_values.get('S6', float('NaN'))
    i6_final_val = calculated_values.get('I6', float('NaN'))
    perf['calc_vov6_abs'] = math.sqrt(2 * i6_final_val / (KP_PMOS * s6_final_val)) if KP_PMOS != 0 and math.isfinite(s6_final_val) and s6_final_val > 0 and math.isfinite(i6_final_val) and i6_final_val >= 0 else float('NaN')

    # Vov7 calculation
    s7_val = calculated_values.get('S7', 0)
    perf['calc_vov7'] = math.sqrt(2 * i6_final_val / (KN_NMOS * s7_val)) if KN_NMOS != 0 and math.isfinite(s7_val) and s7_val > 0 and math.isfinite(i6_final_val) and i6_final_val >= 0 else float('NaN')

    # Warnings
    power_budget = specs.get('Power_max', 10e-3)
    if calc_power > power_budget: warnings.append(f"Calculated Power ({format_value(calc_power*1000, 'mW')}) exceeds budget ({format_value(power_budget*1000, 'mW')}).")
    if icmr_min_calc < spec_icmr_min or icmr_max_calc > spec_icmr_max: warnings.append(f"Std. Calc ICMR [{format_value(icmr_min_calc,'V')}, {format_value(icmr_max_calc,'V')}] doesn't meet spec [{spec_icmr_min}V, {spec_icmr_max}V].")
    target_av_db = specs.get('spec_av_db', 120)
    gain_diff = abs(calc_av_db - target_av_db)
    if not math.isfinite(gain_diff) or gain_diff > 6: warnings.append(f"Calc. Gain ({format_value(calc_av_db,'dB')}) differs >6dB from target ({target_av_db} dB).")
    if z == 0: warnings.append("Z=0 -> λ near zero, calculated gain likely unrealistically high.")
    if not math.isfinite(calculated_values.get('S5', 0)): warnings.append("Step 7 Warning: S5 calculation failed.") # Simplified warning
    if not math.isfinite(calculated_values.get('S6_formula', 0)): warnings.append("Step 10 Warning: S6 (Ratio Formula) calculation failed.")
    results['warnings'] = warnings

    return results
# --- End Replace design_opamp ---
# --- END MODIFIED SECTION ---

# --- GUI Class (OpAmpDesigner) ---
# (No changes needed in GUI class itself for this calculation change)
class OpAmpDesigner(tk.Tk):
    # ... (All GUI methods __init__, initial_setup, _update_spec_display,
    #      _calculate_sr_limit, _handle_param_update, update_design,
    #      on_row_select, show_formula_placeholder remain UNCHANGED
    #      from the previous version that included the typing option
    #      and corrected lambda display) ...
    def __init__(self):
        super().__init__()
        self.title("Interactive Op-Amp Designer")
        self.geometry("950x1180") # Adjusted height for more steps

        self.step_data_map = {}
        self.specs = {}

        content_frame = ttk.Frame(self, padding="10")
        content_frame.pack(fill=tk.BOTH, expand=True)

        # --- 1. Specifications Input Section ---
        spec_input_frame = ttk.LabelFrame(content_frame, text="Specifications Input", padding="10")
        spec_input_frame.pack(fill=tk.X, pady=(0, 5))
        ttk.Label(spec_input_frame, text="Enter Last 3 Digits (XYZ):").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.xyz_var = tk.StringVar()
        self.xyz_entry = ttk.Entry(spec_input_frame, width=10, textvariable=self.xyz_var)
        self.xyz_entry.grid(row=0, column=1, padx=5, pady=5, sticky="w")
        self.calc_button = ttk.Button(spec_input_frame, text="Load Specs & Defaults", command=self.initial_setup)
        self.calc_button.grid(row=0, column=2, padx=10, pady=5)

        # --- 2. Adjustable Parameters Section ---
        slider_frame = ttk.LabelFrame(content_frame, text="Adjustable Parameters", padding="10")
        slider_frame.pack(fill=tk.X, pady=5)
        self.cc_var = tk.DoubleVar(); self.sr_var = tk.DoubleVar(); self.power_budget_var = tk.DoubleVar()
        # Cc Row
        ttk.Label(slider_frame, text="Comp. Cap (Cc):").grid(row=0, column=0, sticky=tk.W, padx=5, pady=2)
        self.cc_slider = ttk.Scale(slider_frame, from_=0.1e-12, to=10e-12, orient=tk.HORIZONTAL, variable=self.cc_var, command=self._handle_param_update)
        self.cc_slider.grid(row=0, column=1, sticky=tk.EW, padx=5, pady=2)
        self.cc_entry = ttk.Entry(slider_frame, width=12, textvariable=self.cc_var)
        self.cc_entry.grid(row=0, column=2, sticky=tk.W, padx=5, pady=2)
        self.cc_label = ttk.Label(slider_frame, text="0.00 pF", width=15); self.cc_label.grid(row=0, column=3, sticky=tk.W, padx=5, pady=2)
        self.cc_entry.bind("<Return>", lambda event: self._handle_param_update(event, source='entry', param_id='cc'))
        self.cc_entry.bind("<FocusOut>", lambda event: self._handle_param_update(event, source='entry', param_id='cc'))
        # SR Row
        ttk.Label(slider_frame, text="Slew Rate (SR):").grid(row=1, column=0, sticky=tk.W, padx=5, pady=2)
        self.sr_slider = ttk.Scale(slider_frame, from_=1e6, to=100e6, orient=tk.HORIZONTAL, variable=self.sr_var, command=self._handle_param_update)
        self.sr_slider.grid(row=1, column=1, sticky=tk.EW, padx=5, pady=2)
        self.sr_entry = ttk.Entry(slider_frame, width=12, textvariable=self.sr_var)
        self.sr_entry.grid(row=1, column=2, sticky=tk.W, padx=5, pady=2)
        self.sr_label = ttk.Label(slider_frame, text="0.0 V/µs", width=15); self.sr_label.grid(row=1, column=3, sticky=tk.W, padx=5, pady=2)
        self.sr_entry.bind("<Return>", lambda event: self._handle_param_update(event, source='entry', param_id='sr'))
        self.sr_entry.bind("<FocusOut>", lambda event: self._handle_param_update(event, source='entry', param_id='sr'))
        # Power Budget Row
        ttk.Label(slider_frame, text="Power Budget:").grid(row=2, column=0, sticky=tk.W, padx=5, pady=2)
        self.power_slider = ttk.Scale(slider_frame, from_=0.5e-3, to=10e-3, orient=tk.HORIZONTAL, variable=self.power_budget_var, command=self._handle_param_update)
        self.power_slider.grid(row=2, column=1, sticky=tk.EW, padx=5, pady=2)
        self.power_entry = ttk.Entry(slider_frame, width=12, textvariable=self.power_budget_var)
        self.power_entry.grid(row=2, column=2, sticky=tk.W, padx=5, pady=2)
        self.power_label = ttk.Label(slider_frame, text="0.0 mW", width=15); self.power_label.grid(row=2, column=3, sticky=tk.W, padx=5, pady=2)
        self.power_entry.bind("<Return>", lambda event: self._handle_param_update(event, source='entry', param_id='power'))
        self.power_entry.bind("<FocusOut>", lambda event: self._handle_param_update(event, source='entry', param_id='power'))

        # --- Start Add Ratio Row ---
        # --- Ratio (I6/I5) Row ---
        ttk.Label(slider_frame, text="I6/I5 Ratio:").grid(row=3, column=0, sticky=tk.W, padx=5, pady=2) # New row 3
        self.ratio_var = tk.DoubleVar(value=5.0) # Default to 5.0
        self.ratio_slider = ttk.Scale(slider_frame, from_=1.0, to=15.0, orient=tk.HORIZONTAL, variable=self.ratio_var, command=self._handle_param_update) # Range 1-15
        self.ratio_slider.grid(row=3, column=1, sticky=tk.EW, padx=5, pady=2)
        self.ratio_entry = ttk.Entry(slider_frame, width=12, textvariable=self.ratio_var)
        self.ratio_entry.grid(row=3, column=2, sticky=tk.W, padx=5, pady=2)
        self.ratio_label = ttk.Label(slider_frame, text="5.00", width=15)
        self.ratio_label.grid(row=3, column=3, sticky=tk.W, padx=5, pady=2)
        self.ratio_entry.bind("<Return>", lambda event: self._handle_param_update(event, source='entry', param_id='ratio'))
        self.ratio_entry.bind("<FocusOut>", lambda event: self._handle_param_update(event, source='entry', param_id='ratio'))
        # --- End Add Ratio Row ---

        # Config
        slider_frame.columnconfigure(1, weight=1)
        # --- Start Add I6 Method Choice ---
        # --- I6 Calculation Method ---
        i6_method_frame = ttk.Frame(slider_frame)
        i6_method_frame.grid(row=4, column=0, columnspan=4, sticky=tk.W, pady=(5,0)) # Spans across columns
        ttk.Label(i6_method_frame, text="I6 Calculation:").pack(side=tk.LEFT, padx=5)
        self.i6_method_var = tk.StringVar(value="ratio") # Default to ratio method
                # --- Start Replace Radiobutton Commands ---
        self.ratio_radio = ttk.Radiobutton(i6_method_frame, text="Use Ratio*I5", variable=self.i6_method_var, value="ratio", command=lambda: self._handle_param_update(event=None, source='radio')) # Pass event=None
        self.ratio_radio.pack(side=tk.LEFT, padx=5)
        self.formula_radio = ttk.Radiobutton(i6_method_frame, text="Use gm6²/ (2 K'p S6_formula)", variable=self.i6_method_var, value="formula", command=lambda: self._handle_param_update(event=None, source='radio')) # Pass event=None
        self.formula_radio.pack(side=tk.LEFT, padx=5)
        # --- End Replace Radiobutton Commands ---
        # --- End Add I6 Method Choice ---
        self.cc_slider.config(state=tk.DISABLED); self.cc_entry.config(state=tk.DISABLED)
        self.sr_slider.config(state=tk.DISABLED); self.sr_entry.config(state=tk.DISABLED)
        self.power_slider.config(state=tk.DISABLED); self.power_entry.config(state=tk.DISABLED)
        # --- Start Disable New Controls --- # ADD THESE LINES
        self.ratio_slider.config(state=tk.DISABLED); self.ratio_entry.config(state=tk.DISABLED)
        self.ratio_radio.config(state=tk.DISABLED); self.formula_radio.config(state=tk.DISABLED)
        # --- End Disable New Controls ---   # END ADDED LINES

        

    # --- 3. Target Specifications Display Section ---
    # (Identical logic, uses spec_list_display)
    # ... (rest of __init__ remains the same) ...
        spec_display_frame = ttk.LabelFrame(content_frame, text="Target Specifications", padding="10")
        spec_display_frame.pack(fill=tk.X, pady=5)
        self.spec_labels = {}
        for i, (text, key, unit) in enumerate(spec_list_display):
            row = 0 + i // 2; col_offset = (i % 2) * 2
            ttk.Label(spec_display_frame, text=text).grid(row=row, column=col_offset, sticky='w', padx=5)
            self.spec_labels[key] = ttk.Label(spec_display_frame, text="-", anchor='w', width=15)
            self.spec_labels[key].grid(row=row, column=col_offset+1, sticky='w', padx=5)

    # --- 4. Design Steps & Formula Section ---
    # (Identical logic, uses latex_formulas)
    # ... (rest of __init__ remains the same) ...
        steps_and_formula_frame = ttk.LabelFrame(content_frame, text="Design Steps & Selected Formula", padding="10")
        steps_and_formula_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        table_frame = ttk.Frame(steps_and_formula_frame)
        table_frame.pack(fill="x", expand=False, pady=(0, 10))
        cols = ("No", "Langkah", "Hasil")
        self.steps_tree = ttk.Treeview(table_frame, columns=cols, show='headings', height=14) # Adjusted height
        for col in cols: self.steps_tree.heading(col, text=col)
        self.steps_tree.column("No", width=50, anchor='center', stretch=False)
        self.steps_tree.column("Langkah", width=350, anchor='w')
        self.steps_tree.column("Hasil", width=200, anchor='e')
        scrollbar_y = ttk.Scrollbar(table_frame, orient="vertical", command=self.steps_tree.yview)
        self.steps_tree.configure(yscrollcommand=scrollbar_y.set)
        self.steps_tree.grid(row=0, column=0, sticky='nsew')
        scrollbar_y.grid(row=0, column=1, sticky='ns')
        table_frame.grid_rowconfigure(0, weight=1); table_frame.grid_columnconfigure(0, weight=1)
        self.steps_tree.bind('<<TreeviewSelect>>', self.on_row_select)
        self.fig = matplotlib.figure.Figure(figsize=(8, 3.5), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=steps_and_formula_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=5, pady=(5, 0))
        self.show_formula_placeholder("Enter XYZ and Load Specs")

    # --- 5. Performance Summary Section ---
    # (Identical logic, uses perf_list_display)
    # ... (rest of __init__ remains the same) ...
        perf_frame = ttk.LabelFrame(content_frame, text="Performance Summary", padding="10")
        perf_frame.pack(fill=tk.X, pady=5)
        self.perf_labels = {}
        for i, (text, key, unit) in enumerate(perf_list_display):
             row = 0 + i // 2; col_offset = (i % 2) * 2
             ttk.Label(perf_frame, text=text).grid(row=row, column=col_offset, sticky='w', padx=5, pady=2)
             # Use the updated key 'calc_vov3_actual' if needed here
             display_key = key if key != 'chosen_vov3_abs' else 'calc_vov3_actual'
             self.perf_labels[display_key] = ttk.Label(perf_frame, text="-", width=15, anchor='w')
             self.perf_labels[display_key].grid(row=row, column=col_offset+1, sticky='w', padx=5, pady=2)


    # --- Status Bar ---
    # (Identical logic)
    # ... (rest of __init__ remains the same) ...
        self.status_label = ttk.Label(content_frame, text="Enter XYZ digits and click 'Load Specs & Defaults'.", relief=tk.SUNKEN, anchor='w', wraplength=900)
        self.status_label.pack(side=tk.BOTTOM, fill='x', pady=(10, 0))

    # --- initial_setup ---
    # (Unchanged from version with typing option)
        # --- Start Replace initial_setup ---
    def initial_setup(self):
        xyz = self.xyz_var.get()
        if len(xyz) != 3 or not xyz.isdigit(): messagebox.showerror("Input Error", "Please enter exactly 3 digits for XYZ."); return
        try:
            x = int(xyz[0]); y = int(xyz[1]); z = int(xyz[2])
            self.specs['CL'] = (10 + x) * 1e-12; self.specs['SR_min'] = x * 1e6; self.specs['VDD'] = 3.0 + y * 0.1; self.specs['Power_max_initial'] = 10e-3
            cc_min = 0.22 * self.specs['CL']; cc_max = 5.0 * self.specs['CL']; cc_default = 0.5 * self.specs['CL']
            self.cc_slider.config(from_=cc_min, to=cc_max, state=tk.NORMAL); self.cc_entry.config(state=tk.NORMAL); self.cc_var.set(cc_default)
            power_min = 0.5e-3; power_max = self.specs['Power_max_initial']; power_default = power_max
            self.power_slider.config(from_=power_min, to=power_max, state=tk.NORMAL); self.power_entry.config(state=tk.NORMAL); self.power_budget_var.set(power_default); self.specs['Power_max'] = power_default
            sr_min = self.specs['SR_min']; initial_sr_max = self._calculate_sr_limit(cc_default, power_default); sr_default = sr_min
            self.sr_slider.config(from_=sr_min, to=initial_sr_max, state=tk.NORMAL); self.sr_entry.config(state=tk.NORMAL); self.sr_var.set(sr_default)

            # --- Enable Ratio and I6 Method Controls ---
            self.ratio_var.set(5.0) # Set default ratio
            self.ratio_slider.config(state=tk.NORMAL); self.ratio_entry.config(state=tk.NORMAL)
            self.i6_method_var.set("ratio") # Default method
            self.ratio_radio.config(state=tk.NORMAL); self.formula_radio.config(state=tk.NORMAL)
            # ---

            self._update_spec_display(xyz)
            self._handle_param_update(event=None) # Trigger initial update
            self.status_label.config(text="Sliders/Entries enabled. Adjust parameters. Press Enter or Tab out of entry box to update.")
        except Exception as e: messagebox.showerror("Setup Error", f"Could not initialize parameters: {e}"); import traceback; traceback.print_exc()
    # --- End Replace initial_setup ---

    # --- _update_spec_display ---
    # (Unchanged from version with corrected lambda display)
    def _update_spec_display(self, xyz):
         if len(xyz) != 3 or not xyz.isdigit(): return
         x = int(xyz[0]); y = int(xyz[1]); z = int(xyz[2])
         lambda_p_value = get_lambda(0.05, z); lambda_n_value = get_lambda(0.01, z)
         lambda_p_str = f"0.05{z}"; lambda_n_str = f"0.01{z}"
         spec_vals = {
             'spec_av_db': 120 + z, 'spec_av_lin': 10**((120+z)/20.0), 'spec_gb_hz': (x*10 + y)*1e6, 'spec_vdd': 3.0 + y*0.1, 'spec_vss': 0.0, 'spec_sr_min': x*1e6,
             'spec_cl': (10 + x)*1e-12, 'spec_lambda_p': lambda_p_str, 'spec_lambda_n': lambda_n_str, 'icmr_spec': f"[-1.0V, 1.5V]"
         }
                  # --- Start Add VSS to self.specs ---
         # Store actual calculated float values for design_opamp
         self.specs['GB'] = spec_vals['spec_gb_hz']
         self.specs['lambda_p'] = lambda_p_value
         self.specs['lambda_n'] = lambda_n_value
         self.specs['icmr_min_spec'] = -1.0
         self.specs['icmr_max_spec'] = 1.5
         self.specs['spec_av_db'] = spec_vals['spec_av_db']
         self.specs['VSS'] = spec_vals['spec_vss'] # ADD THIS LINE
         # VDD, SR_min, CL, Power_max are already in self.specs
         # --- End Add VSS to self.specs ---
         # Update GUI Labels
         for key, label in self.spec_labels.items():
              value = spec_vals.get(key, '-'); unit = next((s[2] for s in spec_list_display if s[1] == key), '')
              if key == 'icmr_spec': label.config(text=value)
              elif key == 'spec_lambda_p' or key == 'spec_lambda_n': label.config(text=f"{value} {unit}") # Display string + unit
              elif isinstance(value, (int, float)):
                   if key == 'spec_sr_min': val_disp, unit_disp = value/1e6, 'V/µs'
                   elif key == 'spec_cl': val_disp, unit_disp = value*1e12, 'pF'
                   elif key == 'spec_gb_hz': val_disp, unit_disp = value/1e6, 'MHz'
                   else: val_disp, unit_disp = value, unit
                   label.config(text=format_value(val_disp, unit_disp))
              else: label.config(text=f"{value}")

    # --- _calculate_sr_limit ---
    # (Unchanged from version with typing option)
        # --- Start Replace _calculate_sr_limit ---
    # --- Helper to Calculate SR Limit --- ## CORRECTED ##
    # (Logic extracted from original on_slider_change, now uses ratio_var)
    def _calculate_sr_limit(self, cc_val, power_val):
        """Calculates the maximum SR based on Cc and Power Budget."""
        vdd = self.specs.get('VDD', 3.3) # Default if not set
        sr_min_req = self.specs.get('SR_min', 1e6) # Default if not set

        # --- Get the ratio from the Tkinter variable ---
        current_ratio = self.ratio_var.get() # Get value from the DoubleVar
        # ---

        sr_max_power = float('inf')
        # Use the current_ratio obtained from the GUI variable
        denominator = cc_val * (1 + current_ratio)
        if vdd > 0 and denominator > 1e-18:
            sr_max_power = power_val / denominator / vdd
        reasonable_sr_abs_max = 1000 * 1e6
        new_sr_limit = max(sr_min_req, min(reasonable_sr_abs_max, sr_max_power))
        return new_sr_limit
    # --- End Replace _calculate_sr_limit ---

    # --- _handle_param_update ---
    # (Unchanged from version with typing option)
        # --- Start Replace _handle_param_update ---
    def _handle_param_update(self, event, source='slider', param_id=None):
        """Handles updates from sliders, entry fields, or radio buttons."""
        if not self.specs: return # Exit if setup not done

        try:
            # 1. Handle Input from Entry Box (if applicable)
            if source == 'entry' and param_id:
                entry_widget, variable, slider_widget = None, None, None
                if param_id == 'cc': entry_widget, variable, slider_widget = self.cc_entry, self.cc_var, self.cc_slider
                elif param_id == 'sr': entry_widget, variable, slider_widget = self.sr_entry, self.sr_var, self.sr_slider
                elif param_id == 'power': entry_widget, variable, slider_widget = self.power_entry, self.power_budget_var, self.power_slider
                elif param_id == 'ratio': entry_widget, variable, slider_widget = self.ratio_entry, self.ratio_var, self.ratio_slider # Handle ratio entry
                else: return

                current_text = entry_widget.get()
                try:
                    value_from_entry = float(current_text)
                    min_val, max_val = slider_widget.cget('from'), slider_widget.cget('to')
                    clamped_value = max(min_val, min(max_val, value_from_entry))
                    variable.set(clamped_value)
                    if abs(clamped_value - value_from_entry) > 1e-12 * abs(value_from_entry):
                         entry_widget.delete(0, tk.END); entry_widget.insert(0, f"{clamped_value:.4e}")
                         self.status_label.config(text=f"{param_id.upper()} value clamped.", foreground="orange")
                except ValueError:
                    messagebox.showwarning("Input Error", f"Invalid number: '{current_text}'. Reverting.")
                    entry_widget.delete(0, tk.END); entry_widget.insert(0, f"{variable.get():.4e}")
                    return

            # 2. Check I6 Calculation Method & Enable/Disable Ratio Controls
            i6_method = self.i6_method_var.get()
            if i6_method == "formula":
                self.ratio_slider.config(state=tk.DISABLED)
                self.ratio_entry.config(state=tk.DISABLED)
            else: # method == "ratio"
                self.ratio_slider.config(state=tk.NORMAL)
                self.ratio_entry.config(state=tk.NORMAL)

            # 3. Common Update Logic
            cc_val = self.cc_var.get()
            sr_val = self.sr_var.get()
            power_val = self.power_budget_var.get()
            ratio_val = self.ratio_var.get() # Get current ratio

            # Store current adjustable parameters in specs dictionary
            self.specs['Power_max'] = power_val
            self.specs['i6_i5_ratio'] = ratio_val # Store ratio
            self.specs['i6_method'] = i6_method # Store chosen method

            # Update display labels
            self.cc_label.config(text=f"{format_value(cc_val*1e12, 'pF')}")
            self.sr_label.config(text=f"{format_value(sr_val/1e6, 'V/µs')}")
            self.power_label.config(text=f"{format_value(power_val*1000, 'mW')}")
            self.ratio_label.config(text=f"{ratio_val:.2f}") # Update ratio label

            # Update SR slider's upper limit
            new_sr_limit = self._calculate_sr_limit(cc_val, power_val) # Doesn't depend on ratio directly here
            current_sr_limit = self.sr_slider.cget('to')
            if abs(current_sr_limit - new_sr_limit) > 1e-9 * abs(new_sr_limit): self.sr_slider.config(to=new_sr_limit)

            # Clamp current SR value
            current_sr = self.sr_var.get(); sr_min_req = self.sr_slider.cget('from')
            final_sr_val = max(sr_min_req, min(new_sr_limit, current_sr))
            if abs(final_sr_val - current_sr) > 1e-12 * abs(current_sr):
                self.sr_var.set(final_sr_val); self.sr_label.config(text=f"{format_value(final_sr_val/1e6, 'V/µs')}")

            # 4. Trigger Design Recalculation
            self.update_design()

        except tk.TclError: pass
        except Exception as e: print(f"Error in _handle_param_update: {e}"); messagebox.showerror("Update Error", f"Error: {e}"); traceback.print_exc()
    # --- End Replace _handle_param_update ---

    # --- on_slider_change (OBSOLETE) ---
    def on_slider_change(self, value_str): self._handle_param_update(event=None)

    # --- update_design ---
    # (Unchanged from version with typing option, but calls modified design_opamp)
    def update_design(self):
        xyz = self.xyz_var.get()
        if not self.specs or len(xyz) != 3 or not xyz.isdigit(): return
        target_cc = self.cc_var.get(); target_sr = self.sr_var.get()
        self.status_label.config(text="Recalculating design...", foreground="blue"); self.update_idletasks()
        try:
            # Pass the GUI's specs dictionary (includes ratio and method)
            results = design_opamp(xyz, target_cc, target_sr, self.specs)
            steps = results['steps']; perf = results['performance']; warnings = results['warnings']
            selected_step_num_str = None; selection = self.steps_tree.selection()
            if selection:
                try: selected_iid = selection[0]; item_values = self.steps_tree.item(selected_iid, 'values'); selected_step_num_str = item_values[0] if item_values and len(item_values) > 0 and item_values[0].isdigit() else None
                except (tk.TclError, IndexError): pass
            self.steps_tree.delete(*self.steps_tree.get_children()); self.step_data_map.clear()
            new_iid_to_select = None
            for step_data in steps:
                step_num = step_data['step']; hasil_formatted = format_value(step_data['Hasil'], step_data['Unit']); langkah = step_data['Langkah']
                if isinstance(step_data['Hasil'], (int, float)): # Formatting
                    if step_data['Unit'] == '' and step_data['Hasil'] > 0 and 'S' in langkah: hasil_formatted = f"{step_data['Hasil']:.2f} (W/L)"
                    elif step_data['Unit'] == 'S': hasil_formatted = format_value(step_data['Hasil'] * 1e3, 'mS', precision=3)
                    elif step_data['Unit'] == 'A': hasil_formatted = format_value(step_data['Hasil'] * 1e6, 'uA', precision=3)
                    elif step_data['Unit'] == 'F': hasil_formatted = format_value(step_data['Hasil'] * 1e12, 'pF', precision=3)
                iid = self.steps_tree.insert('', tk.END, values=(step_num, langkah, hasil_formatted)); self.step_data_map[step_num] = step_data
                if selected_step_num_str and str(step_num) == selected_step_num_str: new_iid_to_select = iid
            if new_iid_to_select: self.steps_tree.selection_set(new_iid_to_select); self.steps_tree.focus(new_iid_to_select); self.steps_tree.see(new_iid_to_select)
            elif self.steps_tree.get_children(): first = self.steps_tree.get_children()[0]; self.steps_tree.selection_set(first); self.steps_tree.focus(first)
            # --- Update performance labels using updated perf list --- ## MODIFIED SECTION ##

            for key, label in self.perf_labels.items():
                value = perf.get(key, '-') # key might be 'calc_vovX_actual' now
                # Find original key and unit from display list for formatting purposes
                original_key = key
                if key == 'calc_vov3_actual': original_key = 'chosen_vov3_abs' # Map back Vov3 key
                if key == 'calc_vov5_actual': original_key = 'chosen_vov5' # Map back Vov5 key

                unit = next((p[2] for p in perf_list_display if p[1] == original_key), '')

                if key == 'calc_power': val_disp, unit_disp = value*1000, 'mW'
                elif 'vov' in key.lower() and isinstance(value, (int, float)): val_disp, unit_disp = abs(value)*1000, 'mV' # Keep Vov in mV
                else: val_disp, unit_disp = value, unit

                if isinstance(val_disp, (int, float)): label.config(text=format_value(val_disp, unit_disp, precision=3))
                else: label.config(text=f"{value}") # Handle NaN, Inf etc.

            # --- END MODIFIED SECTION ---
            final_status = "Design updated."; status_color = "green" # Update Status
            sr_max_allowed = self.sr_slider.cget('to'); sr_min_req = self.sr_slider.cget('from')
            if abs(target_sr - sr_max_allowed) < 1e-9 * target_sr and target_sr > sr_min_req + 1e-9:
                 if not any("SR is limited by Power" in w for w in warnings): warnings.insert(0, f"SR is limited by Power Budget & Cc ({format_value(sr_max_allowed/1e6,'V/us')} max).")
            if warnings: final_status += " Warnings:\n- " + "\n- ".join(warnings); status_color = "orange"
            final_status += "\nSelect row in table to render formula."
            self.status_label.config(text=final_status, foreground=status_color)
            self.on_row_select(None)
        except Exception as e: messagebox.showerror("Calculation Error", f"An error occurred: {e}"); self.status_label.config(text=f"Calculation Error: {e}", foreground="red"); traceback.print_exc(); self.show_formula_placeholder("Calculation Error")

    # --- on_row_select ---
    # (Unchanged from version with typing option)
    def on_row_select(self, event):
        try:
            selection = self.steps_tree.selection()
            if not selection and self.steps_tree.get_children(): selection = (self.steps_tree.get_children()[0],); self.steps_tree.selection_set(selection[0]); self.steps_tree.focus(selection[0])
            if not selection: self.show_formula_placeholder("No steps calculated."); return

            selected_item = selection[0]
            item_values = self.steps_tree.item(selected_item, 'values')
            if item_values and len(item_values) > 0 and item_values[0].isdigit():
                step_number = int(item_values[0])
                step_data = self.step_data_map.get(step_number)
                if not step_data: self.show_formula_placeholder(f"Data not found for step {step_number}"); return

                # Update LaTeX formula for Step 3 to match image if needed
                if step_number == 3:
                    latex_str = r"$\frac{W_3}{L_3}=\frac{W_4}{L_4}=\frac{I_5}{K'_p[V_{DD}-V_{in(max)}-|V_{T3}|+V_{T1}]^2}$"
                else:
                    latex_str = latex_formulas.get(step_number, "Formula Missing")

                variables = step_data.get('vars', {})
                var_display_list = []
                # Includes gm4 now
                                # --- Start Replace var_formats in on_row_select ---

                var_formats = {
                    'CL': ('pF', 1e12), 'Cc': ('pF', 1e12), 'SR': ('V/µs', 1/1e6), 'I1': ('µA', 1e6), 'I3': ('µA', 1e6), 'I4': ('µA', 1e6), 'I5': ('µA', 1e6), 'I6': ('µA', 1e6),
                    'GB': ('MHz', 1/1e6), 'gm1': ('mS', 1e3), 'gm2': ('mS', 1e3), 'gm3': ('mS', 1e3), 'gm4': ('mS', 1e3), 'gm6': ('mS', 1e3),
                    'K\'p': ('µA/V²', 1e6), 'K\'n': ('µA/V²', 1e6), 'beta1': ('A/V²', 1), 'VDD': ('V', 1), 'VSS': ('V', 1),
                    'Vin(max)': ('V', 1), 'Vin(min)': ('V', 1),
                    '|VT3|': ('V', 1), 'VT1': ('V', 1), '|VT03|(max)': ('V', 1), 'VT1(max)': ('V', 1), 'VDS5': ('V', 1),
                    'Vgs1(max)': ('V',1), '|Vov3|': ('mV', 1e3), 'Vov5': ('mV', 1e3),
                    'S1': ('', 1), 'S2': ('', 1), 'S3': ('', 1), 'S4': ('', 1), 'S5': ('', 1),
                    'S6_formula': ('', 1), # Add S6_formula for display if formula method used
                    'S6': ('', 1), 'S7': ('', 1), 'Ratio': ('', 1),
                    'lambda_p': ('V⁻¹', 1), 'lambda_n': ('V⁻¹', 1), 'lambda2': ('V⁻¹', 1), 'lambda4': ('V⁻¹', 1), 'lambda6': ('V⁻¹', 1), 'lambda7': ('V⁻¹', 1),
                }

                # --- End Replace var_formats in on_row_select ---
                xyz = self.xyz_var.get(); z_digit = xyz[2] if len(xyz) == 3 and xyz.isdigit() else '?'
                for name, value in variables.items():
                    unit, scale = var_formats.get(name, ('', 1)); precision = 2 if unit == '' else 3
                    value_fmt = ""
                    if name.startswith('lambda'): # Show string + actual value for lambda
                         unit = 'V⁻¹'; base_lambda_val = 0
                         if name in ['lambda4', 'lambda6', 'lambda_p']: base_str = f"0.05{z_digit}"
                         elif name in ['lambda2', 'lambda7', 'lambda_n']: base_str = f"0.01{z_digit}"
                         else: base_str = "?"
                         value_fmt = f"{base_str} ({format_value(value, '', precision=2)}) {unit}"
                    elif isinstance(value, (int, float)): value_fmt = format_value(value * scale, unit, precision=precision)
                    else: value_fmt = f"{value}{unit}"
                    var_display_list.append(f"{name} = {value_fmt}")
                var_display_string = "  |  ".join(var_display_list)
                # Render
                self.ax.clear(); self.ax.text(0.05, 0.8, latex_str, fontsize=12, va='top', ha='left'); self.ax.text(0.05, 0.4, "Values Used: " + var_display_string, fontsize=9, va='top', ha='left', wrap=True)
                self.ax.get_xaxis().set_visible(False); self.ax.get_yaxis().set_visible(False); self.ax.spines['top'].set_visible(False); self.ax.spines['right'].set_visible(False); self.ax.spines['bottom'].set_visible(False); self.ax.spines['left'].set_visible(False)
                self.fig.tight_layout(pad=0.5); self.canvas.draw()
            else: self.show_formula_placeholder('Error retrieving row data')
        except IndexError: self.show_formula_placeholder('Select a row to view the formula')
        except Exception as e: print(f"Error rendering formula/vars: {e}"); traceback.print_exc(); self.show_formula_placeholder(f'Error rendering:\n{e}', color='red')

    # --- show_formula_placeholder ---
    # (Unchanged from version with typing option)
    def show_formula_placeholder(self, message, color='black'):
        self.ax.clear(); self.ax.text(0.5, 0.5, message, color=color, ha='center', va='center', wrap=True, fontsize=9)
        self.ax.get_xaxis().set_visible(False); self.ax.get_yaxis().set_visible(False); self.ax.spines['top'].set_visible(False); self.ax.spines['right'].set_visible(False); self.ax.spines['bottom'].set_visible(False); self.ax.spines['left'].set_visible(False)
        self.fig.tight_layout(pad=0.5); self.canvas.draw()


# --- Global lists needed by GUI methods --- ## MODIFIED SECTION ##
# (Corrected quotes AND updated Vov3 description)
spec_list_display = [
    ("Target Av (dB):", 'spec_av_db', 'dB'), ("Target Av (Lin):", 'spec_av_lin', ''),
    ("GBW (Hz):", 'spec_gb_hz', 'MHz'),     ("VDD (V):", 'spec_vdd', 'V'),
    ("VSS (V):", 'spec_vss', 'V'),         ("Min SR (V/s):", 'spec_sr_min', 'V/µs'),
    ("CL (F):", 'spec_cl', 'pF'),           ("Lambda_p (1/V):", 'spec_lambda_p', 'V⁻¹'),
    ("Lambda_n (1/V):", 'spec_lambda_n', 'V⁻¹'), ("ICMR Spec (V):", 'icmr_spec', 'V')
]
# --- Start Replace perf_list_display list ---
perf_list_display = [
    ("Calc. Av (Lin):", 'calc_av_lin', ''),   ("Calc. Av (dB):", 'calc_av_db', 'dB'),
    ("Calc. Power (W):", 'calc_power', 'mW'),
    ("Calc. ICMR Min (V):", 'icmr_min_calc', 'V'),("Calc. ICMR Max (V):", 'icmr_max_calc', 'V'),
    ("Vov1 (Calc, V):", 'calc_vov1', 'V'),
    ("|Vov3| (Actual Calc, V):", 'calc_vov3_actual', 'V'), # Updated key and description
    ("|Vov5| (Actual Calc, V):", 'calc_vov5_actual', 'V'), # Updated key and description
    ("|Vov6| (Std Calc, V):", 'calc_vov6_abs', 'V'),
    ("Vov7 (Std Calc, V):", 'calc_vov7', 'V'),
]
# --- End Replace perf_list_display list ---
# --- END MODIFIED SECTION ---


# --- Run the App ---
if __name__ == "__main__":
     app = OpAmpDesigner()
     app.mainloop()
# --- End of code ---