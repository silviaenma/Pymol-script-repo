from pymol import cmd, stored
import tkinter as tk
import math
from tkinter import filedialog, messagebox
from tkinter.scrolledtext import ScrolledText
import os
import csv

pdb_file = None
segmentos_ppii_global = []  # Variable global para guardar segmentos PPII detectados

 # global para controlar si hay proteína cargada

def seleccionar_archivo():
    global pdb_file
    pdb_file = filedialog.askopenfilename(
        title="Selecciona un archivo PDB",
        filetypes=[("Archivos PDB", "*.pdb")]
    )
    if pdb_file:
        cmd.reinitialize()
        cmd.load(pdb_file, "proteina")
        cmd.hide("everything", "proteina")  # Ocultar todo lo visible
        cmd.show("licorice", "proteina")   # Mostrar como licorice
        messagebox.showinfo("Archivo cargado", f"Se cargó:\n{pdb_file}")

def descargar_molecula():
    global pdb_file
    def fetch_pdb():
        global pdb_file
        pdb_id = entry.get().strip()
        if pdb_id:
            try:
                cmd.reinitialize()
                cmd.fetch(pdb_id, name="proteina")
                cmd.hide("everything", "proteina")
                cmd.show("licorice", "proteina")
                # Aquí asignamos pdb_file con un valor para simular archivo cargado
                pdb_file = pdb_id  # Solo guardamos el ID para que funcione la lógica de “selección”
                messagebox.showinfo("Descarga completa", f"Molécula {pdb_id.upper()} descargada y cargada.")
                fetch_window.destroy()
            except Exception as e:
                messagebox.showerror("Error", f"No se pudo descargar {pdb_id}: {e}")
        else:
            messagebox.showwarning("Advertencia", "Por favor ingresa un ID de PDB.")
    
    fetch_window = tk.Toplevel()
    fetch_window.title("Descargar proteína")
    tk.Label(fetch_window, text="ID de la proteína (PDB):").pack(pady=5)
    entry = tk.Entry(fetch_window, width=20)
    entry.pack(pady=5)
    tk.Button(fetch_window, text="Descargar", command=fetch_pdb).pack(pady=5)




def anadir_hidrogenos():
    if not pdb_file:
        messagebox.showwarning("Advertencia", "Primero selecciona un archivo.")
        return
    cmd.h_add("all")
    cmd.sort("all extend 1")
    cmd.show("licorice", "all")
    messagebox.showinfo("Hidrógenos", "Hidrógenos añadidos y mostrados como licorice.")


def eliminar_solventes():
    cmd.remove("solvent")
    messagebox.showinfo("Solventes", "Solventes eliminados.")


def ocultar_side_chains():
    cmd.hide("everything", "proteina and not name N+CA+C+O")
    messagebox.showinfo("Backbone", "Cadenas laterales ocultas (solo backbone).")


def separar_cadenas():
    stored.chains = []
    cmd.iterate("proteina", "stored.chains.append(chain)")
    for cadena in set(stored.chains):
        nuevo_objeto = f"cadena_{cadena}"
        cmd.create(nuevo_objeto, f"proteina and chain {cadena}")


def obtener_angulos_phi_psi_por_cadena(objeto="proteina"):
    resultados = {}
    chains = cmd.get_chains(objeto)
    for chain in chains:
        sel_ca = f"{objeto} and chain {chain} and name CA"
        phipsi = cmd.get_phipsi(sel_ca)
        if not phipsi:
            continue
        for (obj, idx), (phi, psi) in sorted(phipsi.items()):
            if phi is None or psi is None:
                continue
            stored.info = []
            cmd.iterate(f"({obj}`{idx})", "stored.info.append((chain, resn, resi))", space={'stored': stored})
            if not stored.info:
                continue
            ch, resn, resi = stored.info[0]
            resultados[(ch, resi)] = (resn, phi, psi)
    return resultados


def guardar_csv_angulos_phi_psi():
    if not pdb_file:
        messagebox.showwarning("Advertencia", "Primero selecciona un archivo.")
        return

    phi_map = obtener_angulos_phi_psi_por_cadena("proteina")

    stored.res_list = []
    cmd.iterate("proteina and name CA", "stored.res_list.append((chain, resn, resi))")

    datos_csv = [("Cadena", "Residuo", "Número", "Phi", "Psi")]

    for chain, resn, resi in sorted(stored.res_list, key=lambda x: (x[0], int(x[2]))):
        key = (chain, resi)
        if key in phi_map:
            resn_val, phi, psi = phi_map[key]
            datos_csv.append((chain, resn, resi, f"{phi:.2f}", f"{psi:.2f}"))

    # Guardar en CSV con separador ;
    if len(datos_csv) > 1:
        ruta_csv = os.path.join(os.getcwd(), "angulos_phi_psi.csv")
        with open(ruta_csv, mode="w", newline="", encoding="utf-8") as file:
            writer = csv.writer(file, delimiter=";")
            writer.writerows(datos_csv)


def localizar_atomicos_clave_segmentos(segmentos):
    cmd.delete("esfera_*")
    objetos_por_segmento = []
    cmd.h_add()

    for idx, seg in enumerate(segmentos, start=1):
        objetos_segmento = []

        for (resn, resi, chain, _, _) in seg:
            sele_base = f"proteina and chain {chain} and resi {resi}"
            stored.coords = []

            # Átomo CA (carbono alfa)
            stored.ca_coords = []
            cmd.iterate_state(
                1,
                f"{sele_base} and name CA",
                "stored.ca_coords.append((x, y, z))",
                space={'stored': stored}
            )

            # Oxígeno carbonilo
            stored.o_coords = []
            cmd.iterate_state(
                1,
                f"{sele_base} and name O",
                "stored.o_coords.append((x, y, z))",
                space={'stored': stored}
            )

            # Hidrógeno unido al CA (usar nombre fijo "H" en lugar del nombre real)
            stored.h_coords = []
            cmd.iterate_state(
                1,
                f"(neighbor ({sele_base} and name CA)) and elem H",
                "stored.h_coords.append((x, y, z))",  # Solo coordenadas, no nombre
                space={'stored': stored}
            )

            # Crear pseudoatomos con nombres consistentes
            if stored.ca_coords:
                x, y, z = stored.ca_coords[0]
                esfera_name = f"esfera_s{idx}_CA_{resn}_{resi}_{chain}"
                cmd.pseudoatom(esfera_name, pos=[x, y, z])
                cmd.set("sphere_scale", 0.3, esfera_name)
                cmd.color("blue", esfera_name)
                objetos_segmento.append(esfera_name)

            if stored.o_coords:
                x, y, z = stored.o_coords[0]
                esfera_name = f"esfera_s{idx}_O_{resn}_{resi}_{chain}"
                cmd.pseudoatom(esfera_name, pos=[x, y, z])
                cmd.set("sphere_scale", 0.3, esfera_name)
                cmd.color("red", esfera_name)
                objetos_segmento.append(esfera_name)

            if stored.h_coords:
                x, y, z = stored.h_coords[0]
                # NOMBRE FIJADO COMO "H" (en lugar del nombre real del átomo)
                esfera_name = f"esfera_s{idx}_H_{resn}_{resi}_{chain}"
                cmd.pseudoatom(esfera_name, pos=[x, y, z])
                cmd.set("sphere_scale", 0.3, esfera_name)
                cmd.color("green", esfera_name)
                objetos_segmento.append(esfera_name)

        objetos_por_segmento.append(objetos_segmento)

    cmd.group("atomos_clave", "esfera_*")
    return objetos_por_segmento




def calcular_distancias_colindantes(objetos_por_segmento, max_dist=5.0, archivo_salida="distancias_colindantes.txt"):
    tripletas_candidatas = []
    
    # Precalcular coordenadas de todos los objetos
    coordenadas = {}
    for segmento in objetos_por_segmento:
        for obj in segmento:
            coords = cmd.get_coords(obj)
            if coords is not None and len(coords) > 0:
                coordenadas[obj] = coords[0]  # Tomamos la primera posición

    # Funciones de utilidad
    def es_CA(obj_name): return "_CA_" in obj_name
    def es_O(obj_name): return "_O_" in obj_name
    def es_H(obj_name): return any(x in obj_name for x in ["_HA_", "_H_H_", "_H_CA_", "_H_"])

    def extraer_resi_chain(obj_name):
        partes = obj_name.rsplit('_', 2)
        if len(partes) >= 3:
            return partes[1], partes[2]  # resi, chain
        return None, None

    def distancia(p1, p2):
        return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))

    max_dist_sq = max_dist ** 2

    for i in range(len(objetos_por_segmento) - 1):
        seg_actual = objetos_por_segmento[i]
        seg_siguiente = objetos_por_segmento[i + 1]

        ca_actual = [obj for obj in seg_actual if es_CA(obj)]
        o_actual = [obj for obj in seg_actual if es_O(obj)]
        
        ca_siguiente = [obj for obj in seg_siguiente if es_CA(obj)]
        o_siguiente = [obj for obj in seg_siguiente if es_O(obj)]

        mapa_h_actual = {}
        mapa_h_siguiente = {}
        
        for obj in seg_actual:
            if es_H(obj):
                resi, chain = extraer_resi_chain(obj)
                if resi and chain and obj in coordenadas:
                    mapa_h_actual[(resi, chain)] = obj
        
        for obj in seg_siguiente:
            if es_H(obj):
                resi, chain = extraer_resi_chain(obj)
                if resi and chain and obj in coordenadas:
                    mapa_h_siguiente[(resi, chain)] = obj

        for ca in ca_actual:
            if ca not in coordenadas:
                continue
                
            resi_ca, chain_ca = extraer_resi_chain(ca)
            h_asociado = mapa_h_actual.get((resi_ca, chain_ca), None)
            coord_ca = coordenadas[ca]
            
            for o in o_siguiente:
                if o not in coordenadas:
                    continue
                    
                coord_o = coordenadas[o]
                dist_sq = sum((a - b) ** 2 for a, b in zip(coord_ca, coord_o))
                
                if dist_sq < max_dist_sq:
                    dist = math.sqrt(dist_sq)
                    tripletas_candidatas.append((ca, o, h_asociado, dist))

        for o in o_actual:
            if o not in coordenadas:
                continue
                
            coord_o = coordenadas[o]
            
            for ca in ca_siguiente:
                if ca not in coordenadas:
                    continue
                    
                resi_ca, chain_ca = extraer_resi_chain(ca)
                h_asociado = mapa_h_siguiente.get((resi_ca, chain_ca), None)
                coord_ca = coordenadas[ca]
                
                dist_sq = sum((a - b) ** 2 for a, b in zip(coord_o, coord_ca))
                
                if dist_sq < max_dist_sq:
                    dist = math.sqrt(dist_sq)
                    tripletas_candidatas.append((o, ca, h_asociado, dist))

    # Guardar reporte en archivo
    with open(archivo_salida, "w") as f:
        if tripletas_candidatas:
            f.write(f"Tripletas colindantes con distancia < {max_dist:.1f} Å:\n")
            for a1, a2, a3, dist in tripletas_candidatas:
                f.write(f"  {a1} - {a2} - {a3 if a3 else 'Sin H'} : {dist:.2f} Å\n")
        else:
            f.write(f"No se encontraron pares con distancia < {max_dist:.1f} Å.\n")

    return tripletas_candidatas



def calcular_angulos_ca_h_o(tripletas_candidatas, archivo_salida="angulos_ca_h_o.txt"):
    """
    Para cada tripleta (a1, a2, h, distancia) en tripletas_candidatas,
    calcula el ángulo CA–H–O y lo guarda en un .txt **solo** si el ángulo está entre 110° y 180°.

    Parámetros:
     - tripletas_candidatas: lista de tuplas (obj1, obj2, objH, distancia)
     - archivo_salida: ruta del archivo donde se guardan los ángulos

    Retorna:
     - lista de tuplas (CA, H, O, angulo_en_grados) filtrada con ángulos [110, 180]
    """
    angulos = []

    with open(archivo_salida, "w") as f:
        f.write("CA - H - O : Ángulo (grados)\n")
        for obj1, obj2, objH, dist in tripletas_candidatas:
            # Saltar si no hay hidrógeno asociado
            if objH is None:
                continue

            # Determinar cuál es CA y cuál es O
            if "_CA_" in obj1:
                ca, o = obj1, obj2
            else:
                ca, o = obj2, obj1

            # Calcular ángulo CA–H–O
            try:
                ang = cmd.get_angle(ca, objH, o)
            except Exception as e:
                # Si falla PyMOL en la selección, lo informamos y seguimos
                print(f"Error calculando ángulo para {ca}, {objH}, {o}: {e}")
                continue

            # Filtrar solo ángulos entre 110° y 180°
            if 110.0 <= ang <= 180.0:
                angulos.append((ca, objH, o, ang))
                f.write(f"{ca} - {objH} - {o} : {ang:.2f}°\n")

    return angulos


def calcular_y_visualizar_distancias(max_dist=5.0):
    global segmentos_ppii_global
    if not segmentos_ppii_global:
        messagebox.showwarning("Advertencia", "Primero detecta los segmentos PPII.")
        return

    objetos = localizar_atomicos_clave_segmentos(segmentos_ppii_global)
    pares = calcular_distancias_colindantes(objetos, max_dist=max_dist)
    visualizar_distancias_pares(pares)
    calcular_angulos_ca_h_o(pares)



def visualizar_distancias_pares(pares_candidatos):
    """
    Recibe lista de tuplas (atomo1, atomo2, atomo_H, distancia) y crea objetos distancia en PyMOL.
    """
    if not pares_candidatos:
        messagebox.showinfo("Visualización", "No hay pares de átomos para visualizar.")
        return

    cmd.delete("distancia_ppii")

    for i, (at1, at2, _, dist) in enumerate(pares_candidatos, start=1):
        nombre_dist = f"distancia_ppii_{i}"
        cmd.distance(nombre_dist, at1, at2)
        cmd.set("dash_width", 4, nombre_dist)
        cmd.set("dash_length", 0.5, nombre_dist)
        cmd.color("cyan", nombre_dist)

    messagebox.showinfo("Visualización", f"Visualizados {len(pares_candidatos)} pares de distancias colindantes.")



def detectar_segmentos_ppii(objeto="proteina", min_length=3, tol=20.0):
    global segmentos_ppii_global

    if not pdb_file:
        messagebox.showwarning("Advertencia", "Primero selecciona un archivo.")
        return

    phi_psi_map = obtener_angulos_phi_psi_por_cadena(objeto)
    lista_residuos = []
    for (chain, resi), (resn, phi, psi) in phi_psi_map.items():
        try:
            resi_num = int(resi)
        except:
            continue
        lista_residuos.append((chain, resi_num, resn, phi, psi))
    lista_residuos.sort(key=lambda x: (x[0], x[1]))

    segmentos = []
    segmento_actual = []

    def en_rango_ppii(phi, psi, tol=tol):
        return abs(phi + 75) <= tol and abs(psi - 145) <= tol

    for i, (chain, resi, resn, phi, psi) in enumerate(lista_residuos):
        if en_rango_ppii(phi, psi):
            if not segmento_actual:
                segmento_actual.append((resn, resi, chain, phi, psi))
            else:
                _, last_resi, last_chain, _, _ = segmento_actual[-1]
                if chain == last_chain and resi == last_resi + 1:
                    segmento_actual.append((resn, resi, chain, phi, psi))
                else:
                    if len(segmento_actual) >= min_length:
                        segmentos.append(segmento_actual)
                    segmento_actual = [(resn, resi, chain, phi, psi)]
        else:
            if len(segmento_actual) >= min_length:
                segmentos.append(segmento_actual)
            segmento_actual = []
    if len(segmento_actual) >= min_length:
        segmentos.append(segmento_actual)

    if not segmentos:
        messagebox.showinfo("Resultado", "No se encontraron segmentos PPII.")
        return

    cmd.delete("ppii_segmento*")

    salida = "Segmentos candidatos a hélices PPII:\n"
    for idx, seg in enumerate(segmentos, start=1):
        start_resi = seg[0][1]
        end_resi = seg[-1][1]
        chain = seg[0][2]
        salida += f"\nSegmento {idx} (Cadena {chain}, residuos {start_resi}-{end_resi}, longitud {len(seg)}):\n"
        for (resn, resi, _, phi, psi) in seg:
            salida += f"  {resn}-{resi}{chain}: (phi={phi:.1f}, psi={psi:.1f})\n"

        sel_str = f"proteina and chain {chain} and resi {start_resi}-{end_resi}"
        obj_name = f"ppii_segmento_{chain}_{start_resi}_{end_resi}"
        cmd.create(obj_name, sel_str)
        cmd.color("red", obj_name)
        cmd.show("cartoon", obj_name)

    ruta_archivo = os.path.join(os.getcwd(), "segmentos_ppii.txt")
    with open(ruta_archivo, "w") as f:
        f.write(salida)

    # Guardar global
    segmentos_ppii_global = segmentos

    messagebox.showinfo("Éxito", f"{len(segmentos)} segmentos PPII detectados, resaltados.\n"
                                 f"Átomos clave también visualizados en PyMOL.")


def guardar_segmentos_ppii_pdb():
    global segmentos_ppii_global
    if not segmentos_ppii_global:
        messagebox.showwarning("Advertencia", "Primero detecta los segmentos PPII.")
        return

    # Para cada segmento almacenado, construimos el nombre de objeto y hacemos cmd.save
    count = 0
    for seg in segmentos_ppii_global:
        start_resi = seg[0][1]
        end_resi = seg[-1][1]
        chain = seg[0][2]
        obj_name = f"ppii_segmento_{chain}_{start_resi}_{end_resi}"
        # Comprobamos que el objeto efectivamente exista en PyMOL
        if cmd.count_atoms(f"{obj_name}") > 0:
            filename = os.path.join(os.getcwd(), f"{obj_name}.pdb")
            try:
                cmd.save(filename, obj_name)
                count += 1
            except Exception as e:
                print(f"Error guardando {obj_name}: {e}")
        else:
            print(f"Objeto {obj_name} no encontrado en la sesión de PyMOL.")

    if count > 0:
        messagebox.showinfo("Éxito", f"Se guardaron {count} archivos PDB de segmentos PPII:\n"
                                     f"{os.getcwd()}")
    else:
        messagebox.showwarning("Atención", "No se guardó ningún segmento (quizá no existan objetos en PyMOL).")


def convertir_a_selecciones_pymol(pares_con_distancias):
    selecciones = []
    for at1, at2, _ in pares_con_distancias:
        sele1 = f"id {cmd.index(at1)[0][1]}"
        sele2 = f"id {cmd.index(at2)[0][1]}"
        selecciones.append((sele1, sele2))
    return selecciones




def lanzar_interfaz():
    root = tk.Tk()
    root.title("Análisis de Proteína en PyMOL")
    tk.Button(root, text="Seleccionar archivo PDB", command=seleccionar_archivo, width=40).pack(pady=10)
    tk.Button(root, text="Descargar proteina", command=descargar_molecula, width=40).pack(pady=10)
    tk.Button(root, text="Añadir hidrógenos", command=anadir_hidrogenos, width=40).pack(pady=10)
    tk.Button(root, text="Eliminar solventes", command=eliminar_solventes, width=40).pack(pady=10)
    tk.Button(root, text="Ocultar cadenas laterales", command=ocultar_side_chains, width=40).pack(pady=10)
    tk.Button(root, text="Guardar ángulos phi/psi en archivo", command=guardar_csv_angulos_phi_psi, width=40).pack(pady=10)
    tk.Button(root, text="Detectar segmentos PPII y resaltarlos", command=detectar_segmentos_ppii, width=40).pack(pady=10)
    tk.Button(root, text="Guardar segmentos PPII en PDB", command=guardar_segmentos_ppii_pdb, width=40).pack(pady=10)
    tk.Button(root, text="Calcular y visualizar distancias entre CA-O colindantes",command=calcular_y_visualizar_distancias, width=40).pack(pady=10)


    root.mainloop()


lanzar_interfaz()
