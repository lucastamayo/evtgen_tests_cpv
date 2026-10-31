filename = "events_SSDCP2_pythia.txt"

n_signal_B = 0
b0_to_b0 = 0
b0_to_b0bar = 0
b0bar_to_b0bar = 0
b0bar_to_b0 = 0

examples_printed = 0


def process_event(event, iev):
    global n_signal_B, b0_to_b0, b0_to_b0bar
    global b0bar_to_b0bar, b0bar_to_b0, examples_printed

    # pré-calcula pids (MUITO mais rápido)
    pid_map = {i: p["pid"] for i, p in event.items()}

    for idx, p in event.items():

        if abs(p["pid"]) != 511:
            continue

        d1, d2 = p["daughter1"], p["daughter2"]
        if d1 <= 0 or d2 <= 0 or d2 < d1:
            continue

        has_jpsi = False
        has_ks = False

        # loop direto (sem lista!)
        for i in range(d1, d2 + 1):
            pid = pid_map.get(i)
            if pid is None:
                continue
            if abs(pid) == 443:
                has_jpsi = True
            if pid == 310:
                has_ks = True

        if not (has_jpsi and has_ks):
            continue

        # ===== B sinal =====
        n_signal_B += 1

        final_pid = p["pid"]

        # subir linhagem (rápido)
        current = p
        while True:
            m = current["mother1"]
            if m <= 0 or m not in event:
                break
            mother = event[m]
            if abs(mother["pid"]) == 511:
                current = mother
            else:
                break

        initial_pid = current["pid"]

        if initial_pid == 511 and final_pid == 511:
            b0_to_b0 += 1
        elif initial_pid == 511 and final_pid == -511:
            b0_to_b0bar += 1
        elif initial_pid == -511 and final_pid == -511:
            b0bar_to_b0bar += 1
        elif initial_pid == -511 and final_pid == 511:
            b0bar_to_b0 += 1

        # debug leve
        if examples_printed < 5:
            print(f"Event {iev}: final={final_pid}, initial={initial_pid}")
            examples_printed += 1


# ============================
# LOOP STREAMING REAL
# ============================

current_event = {}
iev = -1

with open(filename) as f:
    for line in f:
        line = line.strip()

        if not line:
            continue

        if line.startswith("Event"):
            if current_event:
                process_event(current_event, iev)
            current_event = {}
            iev += 1
            continue

        parts = line.split()
        if len(parts) < 8:
            continue

        try:
            idx = int(parts[0])
            pid = int(parts[1])
            mother1 = int(parts[4])
            daughter1 = int(parts[6])
            daughter2 = int(parts[7])
        except:
            continue

        current_event[idx] = {
            "pid": pid,
            "mother1": mother1,
            "daughter1": daughter1,
            "daughter2": daughter2,
        }

# último evento
if current_event:
    process_event(current_event, iev)


# ============================
# RESULTADOS
# ============================

print("\n===== RESULTADOS =====\n")
print(f"Total Bs sinal: {n_signal_B}")
print(f"B0 -> B0: {b0_to_b0}")
print(f"B0 -> anti-B0: {b0_to_b0bar}")
print(f"anti-B0 -> anti-B0: {b0bar_to_b0bar}")
print(f"anti-B0 -> B0: {b0bar_to_b0}")

if n_signal_B > 0:
    print("\nFrações:")
    print(f"B0 -> B0: {b0_to_b0 / n_signal_B:.4f}")
    print(f"B0 -> anti-B0: {b0_to_b0bar / n_signal_B:.4f}")
    print(f"anti-B0 -> anti-B0: {b0bar_to_b0bar / n_signal_B:.4f}")
    print(f"anti-B0 -> B0: {b0bar_to_b0 / n_signal_B:.4f}")