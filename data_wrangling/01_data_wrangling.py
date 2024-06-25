#pip install pyinaturalist geopandas pandas matplotlib shapely
# Import necessary libraries
import time
from tqdm import tqdm
from pyinaturalist import get_observations
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import matplotlib.pyplot as plt

# Define your search parameters for iNaturalist
params = {
    'place_id': '97394',  # Place ID for the United States
    'taxon_id': '47157',  # Taxon ID for Lepidoptera
    'quality_grade': 'research',  # Research grade observations
    'geo': True,
    'd1': '2020-01-01',  # Date range start (modify as needed)
    'd2': '2024-04-30',  # Date range end (modify as needed)
    'per_page': 1000000000,  # Number of results per page
}

# Fetch observations from iNaturalist
response = get_observations(**params)

# Extract relevant data
results = response['results']
data = []
for result in results:
    if result.get('geojson'):
        data.append({
            'id': result['id'],
            'species': result['taxon']['name'] if result.get('taxon') else None,
            'latitude': result['geojson']['coordinates'][1],
            'longitude': result['geojson']['coordinates'][0],
            'date': result['observed_on']
        })

# Create a DataFrame
inat_df = pd.DataFrame(data)

# Convert to a GeoDataFrame
inat_gdf = gpd.GeoDataFrame(
    inat_df,
    geometry=[Point(xy) for xy in zip(inat_df.longitude, inat_df.latitude)],
    crs='EPSG:4326'
)

# Display the first few rows
print(inat_gdf.head())

# Load the greenspace shapefile (ensure all associated files are in the same directory)
greenspace_gdf = gpd.read_file("/Users/jennycheung/Documents/phd_study/iNat_project/data/USA_Parks/USA Parks.shp")

# Ensure both GeoDataFrames use the same CRS
inat_gdf = inat_gdf.to_crs(greenspace_gdf.crs)

# Re-project to a suitable UTM zone (let's use UTM Zone 14N as an example)
projected_crs = "EPSG:32614"  # UTM Zone 14N
greenspace_gdf = greenspace_gdf.to_crs(projected_crs)

# Time the buffering process with a progress meter
start_time = time.time()
greenspace_gdf['geometry'] = [geom.buffer(100) for geom in tqdm(greenspace_gdf['geometry'])]
end_time = time.time()

# Calculate and print the elapsed time
elapsed_time = end_time - start_time
print(f"Buffering took {elapsed_time:.2f} seconds")

# Re-project back to the original CRS if needed
greenspace_gdf = greenspace_gdf.to_crs(inat_gdf.crs)

# Perform spatial join to find iNaturalist points within the greenspace buffers
inat_within_buffers = gpd.sjoin(inat_gdf, greenspace_gdf, predicate='within')

# Plot the result
#base = greenspace_gdf.plot(color='green', edgecolor='black', alpha=0.5)
#inat_within_buffers.plot(ax=base, color='red', markersize=5)
#plt.title("iNaturalist Observations within Greenspace Buffers")

# Save the plot as an image file
#plt.savefig("/Users/jennycheung/Documents/phd_study/iNat_project/plot/inat_within_buffers_plot.png")

# Show the plot
#plt.show()

# Output the results to various formats
# Output to CSV

start_time = time.time()
inat_within_buffers.to_csv("/Users/jennycheung/Documents/phd_study/iNat_project/output_data/inat_within_buffers.csv", index=False)
csv_elapsed_time = time.time() - start_time
print(f"CSV output took {csv_elapsed_time:.2f} seconds")

# Output to GeoJSON
start_time = time.time()
inat_within_buffers.to_file("/Users/jennycheung/Documents/phd_study/iNat_project/output_data/inat_within_buffers.geojson", driver="GeoJSON")
geojson_elapsed_time = time.time() - start_time
print(f"GeoJSON output took {geojson_elapsed_time:.2f} seconds")

# Output to Shapefile
start_time = time.time()
inat_within_buffers.to_file("/Users/jennycheung/Documents/phd_study/iNat_project/output_data/inat_within_buffers.shp")
shapefile_elapsed_time = time.time() - start_time
print(f"Shapefile output took {shapefile_elapsed_time:.2f} seconds")
