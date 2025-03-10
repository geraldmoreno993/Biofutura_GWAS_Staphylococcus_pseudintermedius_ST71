import os
import subprocess
from glob import glob

os.chdir("/home/ubigem/Documentos/SP_segundo_intento")
# Detecta todos los archivos .fasta en la carpeta actual
genome_files = glob("*.fasta")

# Verifica si hay archivos .fasta
if not genome_files:
    print("No se encontraron archivos .fasta en la carpeta actual.")
    exit()

# Archivo de salida para guardar los resultados
output_file = "resultados.txt"

# Comienza con un encabezado en el archivo de salida
with open(output_file, "w") as out:
    out.write("Genome1\tGenome2\tANI\n")

# Comparar todos los pares de genomas
for i, genome1 in enumerate(genome_files):
    for j, genome2 in enumerate(genome_files):
        if i < j:  # Evita repetir comparaciones
            try:
                # Ejecuta OrthoANI para el par actual
                result = subprocess.run(
                    ["orthoani", "-q", genome1, "-r", genome2],
                    capture_output=True,
                    text=True,
                    check=True,
                )

                # Extrae el valor de ANI del resultado
                ani_value = result.stdout.strip().split()[-1]

                # Escribe los resultados en el archivo
                with open(output_file, "a") as out:
                    out.write(f"{genome1}\t{genome2}\t{ani_value}\n")
                print(f"Comparación completada: {genome1} vs {genome2}, ANI: {ani_value}")

            except subprocess.CalledProcessError as e:
                print(f"Error al procesar {genome1} y {genome2}: {e}")

