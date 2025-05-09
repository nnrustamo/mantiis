import pandas as pd
import matplotlib.pyplot as plt

file_path = 'run_log_results.xlsx'
df = pd.read_excel(file_path)

labels = ['(a)', '(b)', '(c)']

def plot_runtime(threads, label):
    thread_data = df[df['threads'] == threads]
    mg_true = thread_data[thread_data['MultiGrid'] == True]
    mg_false = thread_data[thread_data['MultiGrid'] == False]

    plt.rc('font', family='Times New Roman')

    plt.figure(figsize=(5, 5))
    plt.text(0.05, 0.95, label, transform=plt.gca().transAxes, fontsize=16, fontname="Times New Roman", verticalalignment='top', horizontalalignment='left')
    plt.plot(mg_true['Nx'], mg_true['Time Per Iteration (s)']*1000, label='Multi-block LB', color='black', marker='o', markersize=6)
    plt.plot(mg_false['Nx'], mg_false['Time Per Iteration (s)']*1000, label='Uniform grid LB', color='red', marker='D', markersize=6)
    plt.xlabel('Grid size', fontsize=14)
    plt.ylabel('Time per iteration (ms)', fontsize=14)
    plt.ylim([0, 500])
    plt.legend(fontsize=14)
    plt.grid(False)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # plt.savefig(r"{}.png".format(threads))
    plt.tight_layout()
    # plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.7, color='gray')
    # t.show()


plot_runtime(1, labels[0])
plt.savefig("1.svg")
plot_runtime(4, labels[1])
plt.savefig("4.svg")
plot_runtime(8, labels[2])
plt.savefig("8.svg")
