#To use this plot script, it needs the output coming from detective.py
# Load the data from the CSV file
# Load the data from the CSV file
df = pd.read_csv("overlap_results.csv", sep =",")
# Create a scatter plot
plt.figure(figsize=(10, 6))

# Loop through unique Pathogenicity_Island1 values to create separate scatter points
for key, grp in df.groupby(['Pathogenicity_Island1']):
    plt.scatter(grp['Pathogenicity_Island2'], grp['Overlap_Size'], label=key)

plt.title('Overlap Size between Pathogenicity Islands')
plt.xlabel('Pathogenicity Island 2')
plt.ylabel('Overlap Size')
plt.xticks(rotation=45)
plt.legend(title='Pathogenicity Island 1')
plt.grid()
plt.tight_layout()
plt.savefig("overlap_size_plot.pdf")  # Save as PDF
plt.show()
