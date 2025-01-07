import pandas as pd

# Extraer el identificador de muestra (SRRXXXXXXXX) de la columna isoform
def extract_sample_id(isoform):
    match = re.search(r"SRR\d+", isoform)
    return match.group(0) if match else None

# Rutas archivos de SQANTI3-QC
classification_file = "/home/vgaya/estudios_in_silico/ONT/isoform_sqanti3_illumina_2/isoform_illumina_classification.txt"
input_count = "/home/vgaya/estudios_in_silico/ONT/flair_out_2/align.counts_matrix.tsv.counts.tsv"

# Salida
output_count_matrix = "/home/vgaya/estudios_in_silico/ONT/isoform_counts_matrix.tsv"

# Leer los archivos
df1 = pd.read_csv(input_count, sep="\t")
df2 = pd.read_csv(classification_file, sep="\t", dtype={'isoform': str}, low_memory=False)

# Extraer las columnas necesarias del segundo archivo
df2 = df2[['isoform', 'structural_category', 'associated_gene']]

# Crear un diccionario para asociar cada isoform con sus valores
isoform_dict = df2.groupby('isoform').agg({
    'structural_category': lambda x: ";".join(map(str, x.unique())),
    'associated_gene': lambda x: ";".join(map(str, x.unique()))
}).to_dict('index')

# Inicializar nuevas columnas en el primer DataFrame
df1['structural_category'] = ""
df1['associated_gene'] = ""

# Buscar coincidencias y asignar valores
for i, row in df1.iterrows():
    ids = row['ids']
    matches = [key for key in isoform_dict.keys() if key in ids]
    if matches:
        df1.at[i, 'structural_category'] = ";".join(
                isoform_dict[match]['structural_category'] for match in matches)
        df1.at[i, 'associated_gene'] = ";".join(
                isoform_dict[match]['associated_gene'] for match in matches)

# Guardar el archivo de salida
df1.to_csv(output_count_matrix, sep="\t", index=False)

