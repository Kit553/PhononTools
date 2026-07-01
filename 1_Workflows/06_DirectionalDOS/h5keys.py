#!/usr/bin/env python3

import sys
import h5py

def print_keys(name, obj):
    print("/" + name)

def main():
    if len(sys.argv) != 2:
        print("Usage: python print_h5_keys.py file.h5")
        sys.exit(1)

    with h5py.File(sys.argv[1], "r") as h5:
        print("Top-level keys:")
        for key in h5.keys():
            print("/" + key)

        print("\nAll keys:")
        h5.visititems(print_keys)

if __name__ == "__main__":
    main()
