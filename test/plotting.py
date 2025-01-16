import pandas as pd
import matplotlib.pyplot as plt

file_path = 'run_log_results.xlsx'
df = pd.read_excel(file_path)

def plot_runtime(threads):
    thread_data = df[df['threads'] == threads]
    mg_true = thread_data[thread_data['MultiGrid'] == True]
    mg_false = thread_data[thread_data['MultiGrid'] == False]

    plt.rc('font', family='Times New Roman')

    plt.figure(figsize=(8, 6))
    plt.plot(mg_true['Nx'], mg_true['Time Per Iteration (s)']*1000, label='Multi-block', color='black', marker='o', markersize=6)
    plt.plot(mg_false['Nx'], mg_false['Time Per Iteration (s)']*1000, label='Single-block', color='red', marker='D', markersize=6)
    plt.xlabel('Nx', fontsize=12)
    plt.ylabel('Time Per Iteration (ms)', fontsize=12)
    plt.legend(fontsize=12)
    plt.grid(False)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(r"{}.png".format(threads))
    # t.show()


plot_runtime(1)
plt.savefig("1.png")
plot_runtime(4)
plt.savefig("4.png")
plot_runtime(8)
plt.savefig("8.png")
