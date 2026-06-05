# calculate_average.py
import sys

def calculate_average_third_column(file_path):
    total = 0
    count = 0
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    value = float(parts[2])
                    total += value
                    count += 1
                except ValueError:
                    continue  # Skip lines with non-numeric third column
    return total / count if count > 0 else None

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python calculate_average.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]
    average = calculate_average_third_column(file_path)
    if average is not None:
        print(f"{average}")
    else:
        print("0")  # 或者其他默认值，以防没有足够的数据