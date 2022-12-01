from numba import cuda


cc_cores_per_SM_dict = {
    (2, 0): 32,
    (2, 1): 48,
    (3, 0): 192,
    (3, 5): 192,
    (3, 7): 192,
    (5, 0): 128,
    (5, 2): 128,
    (6, 0): 64,
    (6, 1): 128,
    (7, 0): 64,
    (7, 5): 64,
    (8, 0): 64,
    (8, 6): 128,
    (8, 9): 128,
    (9, 0): 128,
    }


def cuda_cores_status():
    total_cores = 0
    sm_number = 0
    compute_capability = 0
    if cuda.detect():
        device = cuda.get_current_device()
        sm_number = getattr(device, "MULTIPROCESSOR_COUNT")
        compute_capability = device.compute_capability
        cores_per_sm = cc_cores_per_SM_dict.get(compute_capability)
        total_cores = cores_per_sm * sm_number
    return total_cores, sm_number, compute_capability


if __name__ == "__main__":
    cuda_cores_status()
