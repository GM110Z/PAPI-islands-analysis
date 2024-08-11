#divide plasmid sequences by sizes 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# Sample data in a DataFrame (replace this with your actual data)
df = pd.read_csv(sys.argv[1] sep="\t")


# Define the size categories
def categorize_size(size):
    if 1000 <= size <= 5000:
        return 'Very Small Plasmids (1-5 kb)'
    elif 5001 <= size <= 10000:
        return 'Small Plasmids (5-10 kb)'
    elif 10001 <= size <= 20000:
        return 'Medium Plasmids (10-20 kb)'
    elif 20001 <= size <= 50000:
        return 'Medium-Large Plasmids (20-50 kb)'
    elif 50001 <= size <= 100000:
        return 'Large Plasmids (50-100 kb)'
    elif 100001 <= size <= 400000:
        return 'Very Large Plasmids (100-300 kb)'
    else:
        return 'Out of Range'

# Apply the categorization to the DataFrame
df['Size_Category'] = df['Size_bp'].apply(categorize_size)

# Display the DataFrame with the new column
df.to_csv('categorized_plasmid_data.csv', index=False)


# Load your data
df2 = pd.read_csv('categorized_plasmid_data.csv')  # Replace with your file path

# Descriptive statistics
desc_stats = df2.groupby('Size_Category')['Size_bp'].describe()
print(desc_stats)

# Distribution plot
plt.figure(figsize=(10, 6))
sns.histplot(df2['Size_bp'], bins=50, kde=True)
plt.title('Distribution of Plasmid Sizes')
plt.xlabel('Size (bp)')
plt.ylabel('Frequency')
plt.show()
