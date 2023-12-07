import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
import pymysql
import re

db = pymysql.connect(host="localhost", user="curcira", passwd="bio466", db="curcira")
# Create a Cursor object to execute queries.
cur = db.cursor()
cur.execute("SELECT * FROM sequences")
class FastaParser:
    @staticmethod
    def read_fasta(fasta_file_name):
        seq = ''
        seq_name, seq_description, sequence = ([], [], [])
        seq_count = 0
        with open(fasta_file_name) as fh:
            for line in fh:
                if line[0] == ">":
                    def_line = line[1: -1]
                    def_line_split = def_line.split(" ", 1)
                    name = def_line_split[0]
                    seq_name.append(name)
                    desc = def_line_split[1]
                    seq_description.append(desc)
                    if seq != '':
                        sequence.append(seq)
                        seq = ''
                        seq_count += 1
                else:
                    seq += line.replace('\n', '')
            if seq != '':
                sequence.append(seq)
                seq_count += 1

        return seq_name, seq_description, sequence, seq_count

class SequenceAnalyzer:
    @staticmethod
    def gcContent(Seq):
        GC_content = (Seq.count("G") + Seq.count("C")) / len(Seq) * 100
        return round(GC_content, 2)

    def find_cpg_islands(sequence):
        cpg_islands = []
        min_length = 6
        min_gc_content = 50
        min_obs_exp_ratio = 0.6

        # Scan for regions that may be CpG islands
        start = 0
        while start < len(sequence):
            if sequence[start].upper() in ['C', 'G']:
                end = start + 1
                while end < len(sequence) and sequence[end].upper() in ['C', 'G']:
                    end += 1
                island_seq = sequence[start:end]

                gc_content = SequenceAnalyzer.gcContent(island_seq)
                cpg_ratio = SequenceAnalyzer.calculate_cpg_ratio(island_seq)

                if len(island_seq) >= min_length and gc_content >= min_gc_content and cpg_ratio >= min_obs_exp_ratio:
                    cpg_islands.append((start, end))

                start = end
            else:
                start += 1

        return cpg_islands
        print(f"Found {len(cpg_islands)} CpG islands.")

    def calculate_cpg_ratio(fragment):
        c_count = fragment.count('C')
        g_count = fragment.count('G')
        cg_count = fragment.count('CG')
        expected_cg_count = (c_count * g_count) / len(fragment)
        ratio = cg_count / expected_cg_count if expected_cg_count > 0 else 0
        return ratio

    def find_homopolymers(sequence, min_length=5):
        homopolymers = []
        pattern = re.compile(r'(A+|T+|G+|C+)')
        for match in pattern.finditer(sequence):
            if len(match.group()) >= min_length:
                homopolymers.append((match.group(), match.start(), match.end()))

        return homopolymers


class FastaViewer(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        self.title("Fasta Viewer")
        self.geometry("1500x1500")

        self.init_components()
        self.display_existing_data()

    def init_components(self):
        # Choose File Button
        choose_file_button = tk.Button(self, text="Choose File", command=self.choose_file)
        choose_file_button.pack(pady=10)


        self.file_name_label = tk.Label(self, text="No file loaded")  # Create the label
        self.file_name_label.pack(padx=5, pady=5)

        # First Panel (Table)
        self.create_table_panel()

        # Second Panel
        self.name_desc_text = tk.Text(self, wrap=tk.WORD, width=80, height=1)
        self.name_desc_text.pack(padx=5, pady=10)

        # Button Frame
        button_frame = tk.Frame(self)
        button_frame.pack(side=tk.TOP, fill=tk.X, pady=10)

        # Spacer Button inside the Button Frame
        spacer_button = tk.Button(button_frame, text="Display Sequence with Spacer", command=self.spacer)
        spacer_button.pack(side=tk.LEFT, padx=5, pady=5)

        # Motif Search Button inside the Button Frame
        motif_search_button = tk.Button(button_frame, text="Motif Search", command=self.search_motif)
        motif_search_button.pack(side=tk.LEFT, padx=5, pady=5)

        # Entry for Motif inside the Button Frame
        self.motif_entry = tk.Entry(button_frame, width=20)
        self.motif_entry.pack(side=tk.LEFT, padx=5, pady=5)

        # CpG Island Detection Button inside the Button Frame
        cpg_island_button = tk.Button(button_frame, text="Detect CpG Islands", command=self.detect_cpg_islands)
        cpg_island_button.pack(side=tk.LEFT, padx=5, pady=5)

        # Homopolymer Detection Button inside the Button Frame
        homopolymer_button = tk.Button(button_frame, text="Detect Homopolymers", command=self.detect_homopolymers)
        homopolymer_button.pack(side=tk.LEFT, padx=5, pady=5)

        # Clear button
        clear_button = tk.Button(button_frame, text="Clear data", command=self.clear_file)
        clear_button.pack(side=tk.LEFT, padx=5, pady=5)

        # Third Panel (Sequence Display)
        self.sequence_text = tk.Text(self, wrap=tk.WORD, width=85, height=25)
        self.sequence_text.pack(side=tk.LEFT, padx=10, pady=10)

        seq_scrollbar = ttk.Scrollbar(self, orient="vertical", command=self.sequence_text.yview)
        seq_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.sequence_text.configure(yscrollcommand=seq_scrollbar.set)

        self.motif_text = tk.Text(self, wrap=tk.WORD, width=25, height=25)
        self.motif_text.pack(side=tk.LEFT, padx=10, pady=10)

        motif_scrollbar = ttk.Scrollbar(self, orient="vertical", command=self.motif_text.yview)
        motif_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.motif_text.configure(yscrollcommand=motif_scrollbar.set)

        self.hp_text = tk.Text(self, wrap=tk.WORD, width=25, height=25)
        self.hp_text.pack(side=tk.LEFT, padx=10, pady=10)

        hp_scrollbar = ttk.Scrollbar(self, orient="vertical", command=self.hp_text.yview)
        hp_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.hp_text.configure(yscrollcommand=hp_scrollbar.set)

        self.cpg_text = tk.Text(self, wrap=tk.WORD, width=25, height=25)
        self.cpg_text.pack(side=tk.LEFT, padx=10, pady=10)

        cpg_scrollbar = ttk.Scrollbar(self, orient="vertical", command=self.cpg_text.yview)
        cpg_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.cpg_text.configure(yscrollcommand=cpg_scrollbar.set)


    def create_table_panel(self):
        self.tree = ttk.Treeview(self, columns=("Sequence Name", "Description", "Length", "A", "T", "G", "C", "GC"), show="headings")

        self.tree.heading("Sequence Name", text="Sequence Name")
        self.tree.heading("Description", text="Description")
        self.tree.heading("Length", text="Length")
        self.tree.heading("A", text="A count")
        self.tree.heading("T", text="T count")
        self.tree.heading("G", text="G count")
        self.tree.heading("C", text="C count")
        self.tree.heading("GC", text="GC content")

        self.tree.column("Sequence Name", width=200, anchor=tk.W)
        self.tree.column("Description", width=400, anchor=tk.W)  # Adjusted width
        self.tree.column("Length", width=100, anchor=tk.CENTER)
        self.tree.column("A", width=50, anchor=tk.CENTER)
        self.tree.column("T", width=50, anchor=tk.CENTER)
        self.tree.column("G", width=50, anchor=tk.CENTER)
        self.tree.column("C", width=50, anchor=tk.CENTER)
        self.tree.column("GC", width=50, anchor=tk.CENTER)

        self.tree.pack(padx=10, pady=10, fill=tk.BOTH)

        # Vertical Scrollbar
        scrollbar = ttk.Scrollbar(self, orient="vertical", command=self.tree.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.tree.configure(yscrollcommand=scrollbar.set)

        self.tree.bind("<ButtonRelease-1>", self.on_tree_click)

    def display_existing_data(self):
        self.seqDict = {}
        for row in cur.fetchall():
            sequence_name, sequence_description, sequence, sequence_length = row
            A_count = sequence.count("A")
            T_count = sequence.count("T")
            G_count = sequence.count("G")
            C_count = sequence.count("C")
            GC_content = SequenceAnalyzer.gcContent(sequence)
            self.seqDict[sequence_name] = sequence

            self.tree.insert("", "end", values=(sequence_name, sequence_description, sequence_length, A_count, T_count, G_count, C_count, GC_content))
        

    def on_tree_click(self, event):
        self.cpg_text.delete(1.0, tk.END)
        self.hp_text.delete(1.0, tk.END)
        self.motif_text.delete(1.0, tk.END)
        selected_item = self.tree.selection()
        if selected_item:
            sequence_data = self.tree.item(selected_item, 'values')
            self.seq_name = sequence_data[0]
            name_and_desc = self.get_name_desc(self.seq_name, sequence_data[1])
            self.display_name_desc(name_and_desc)
            full_sequence = self.get_full_sequence(sequence_data[0])
            self.display_full_sequence(full_sequence)

    def choose_file(self):
        file_path = filedialog.askopenfilename(title="Choose a Fasta file", filetypes=[("Fasta Files", "*.fasta")])

        if file_path:
            self.display_fasta_content(file_path)

    def display_fasta_content(self, file_path):
        seq_name, seq_description, sequence, seq_count = FastaParser.read_fasta(file_path)

        for i in range(seq_count):
            sequence_name = seq_name[i]
            self.seqDict[sequence_name] = sequence[i]
            sequence_description = seq_description[i]
            sequence_length = len(sequence[i])
            A_count = sequence[i].count("A")
            T_count = sequence[i].count("T")
            G_count = sequence[i].count("G")
            C_count = sequence[i].count("C")
            GC_content = SequenceAnalyzer.gcContent(sequence[i])

            self.tree.insert("", "end", values=(sequence_name, sequence_description, sequence_length, A_count, T_count, G_count, C_count, GC_content))

        self.sequence_text.delete(1.0, tk.END)  # Clear the sequence display
        self.file_name_label["text"] = "File: " + file_path


    def get_name_desc(self, sequence_name, sequence_description):
        formatted_name_desc = f">{sequence_name} {sequence_description}"
        return formatted_name_desc

    def display_name_desc(self, name_desc):
        self.name_desc_text.delete(1.0, tk.END)
        self.name_desc_text.insert(tk.END, name_desc)

    def get_full_sequence(self, sequence_name):
        sequence_to_show = self.seqDict.get(sequence_name)
        formatted_seq = ""
        for i in range(0, len(sequence_to_show), 80):
            formatted_seq += "  " + sequence_to_show[i:i+80] + " \n"

        return formatted_seq

    def display_full_sequence(self, sequence):
        self.sequence_text.delete(1.0, tk.END)
        self.sequence_text.insert(tk.END, sequence)

    def display_motif(self, info):
        self.motif_text.delete(1.0, tk.END)
        self.motif_text.insert(tk.END, info)

    def display_hp(self, info):
        self.hp_text.delete(1.0, tk.END)
        self.hp_text.insert(tk.END, info)

    def display_cpg(self, info):
        self.cpg_text.delete(1.0, tk.END)
        self.cpg_text.insert(tk.END, info)


    def search_motif(self):
        target = self.motif_entry.get().upper()
        sequence = self.seqDict.get(self.seq_name)
        motifs = self.find_motifs(target, sequence)
        motifs_print_ready = self.print_motifs(motifs)
        motif_seq = self.convert_motifs_to_lowercase(motifs)
        self.display_full_sequence(motif_seq)
        self.display_motif(motifs_print_ready)



    def find_motifs(self, target, seq):
        motifs = []
        result = re.finditer(target, seq)

        for i in result:
            motifs.append((i.start(), i.end()))
        return motifs

    def print_motifs(self, motifs):
        motifs_print = "Motifs: \n"
        for i in range(0, len(motifs)):
            motif = motifs[i]
            start = int(motif[0])
            end = int(motif[1])
            length = end - start
            motifs_print += str(i + 1) + "=" + str(start) + "-" + str(end) + "_" + str(length) + "\n"
        return motifs_print

    def convert_motifs_to_lowercase(self, motifs=None):
        if motifs:
            sequence = self.seqDict.get(self.seq_name)
            for i in range(0, len(motifs)):
                motif = motifs[i]
                start = int(motif[0])
                end = int(motif[1])
                lowercase = sequence[start:end].lower()
                if start == 0:
                    sequence = lowercase + sequence[end:]
                else:
                    sequence = sequence[0:start] + lowercase + sequence[end:]
            formatted_seq = ""
            for i in range(0, len(sequence), 80):
                formatted_seq += "  " + sequence[i:i+80] + " \n"
        return formatted_seq

    def spacer(self):
        sequence = self.seqDict.get(self.seq_name)
        formatted_sequence = self.format_sequence_string(sequence)
        self.display_full_sequence(formatted_sequence)

    def format_sequence_string(self, sequence):
        sequence = sequence.replace(" \n", "")
        sequence = sequence.replace("  ", "")

        header1_pt2 = ""
        header1_pt1 = (" " * 14) + "1"
        for a in range(2, (70 // 10) + 1):
            header1_pt2 += " " * 10 + str(a)

        header1 = header1_pt1 + header1_pt2
        header2 = "Line " + "1234567890 " * (70 // 10)

        line_count = 0
        formatted_seq = header1 + "\n" + header2
        for j in range(0, len(sequence), 70):
            seq_padded = ""
            seq_line = sequence[j:j + 70]
            line_count += 1
            formatted_seq += "\n{:4d}".format(line_count)
            for b in range(0, 70, 10):
                seq_padded += " " + seq_line[b:b + 10]
            formatted_seq += seq_padded

        return formatted_seq

    def detect_cpg_islands(self):
        sequence = self.sequence_text.get(1.0, tk.END).strip().upper()
        if sequence:
            print("Analyzing sequence for CpG islands...")
            cpg_islands = SequenceAnalyzer.find_cpg_islands(sequence)
            print(f"Found CpG islands: {cpg_islands}")
            cpg_display_format = "CpG Islands"
            for i in cpg_islands:
                start = int(i[0])
                end = int(i[1])
                length = end - start
                cpg_display_format += "\n" + str(start) + "-" + str(end) + "_" + str(length)

            self.cpg_text.insert(tk.END, cpg_display_format)


            if cpg_islands:
                self.highlight_cpg_islands(sequence, cpg_islands)
            else:
                print("No CpG Islands found.")
                self.sequence_text.insert(tk.END, "\nNo CpG Islands found.")
        else:
            print("No sequence available for CpG Island detection.")
            self.sequence_text.insert(tk.END, "\nNo sequence available for CpG Island detection.")

    def highlight_cpg_islands(self, sequence, cpg_islands):
        full_text = self.sequence_text.get("1.0", tk.END)
        print(f"Full text length: {len(full_text)}")  # Diagnostic print

        # Clear previous highlights
        self.sequence_text.tag_remove("cpg_island", "1.0", tk.END)

        for start, end in cpg_islands:
            island_text = sequence[start:end]
            island_start_index = full_text.find(island_text)
            island_end_index = island_start_index + len(island_text)

            print(
                f"Highlighting from {island_start_index} to {island_end_index} for island text: {island_text}")

            start_line_char = self.convert_flat_index_to_line_char(island_start_index, full_text)
            end_line_char = self.convert_flat_index_to_line_char(island_end_index, full_text)

            # Apply highlighting
            self.sequence_text.tag_add("cpg_island", start_line_char, end_line_char)

        self.sequence_text.tag_config("cpg_island", background="yellow")

    def convert_flat_index_to_line_char(self, flat_index, text):
        lines = text.split('\n')
        current_index = 0
        for line_number, line in enumerate(lines, start=1):
            if current_index + len(line) + 1 > flat_index:
                char_index = flat_index - current_index
                return f"{line_number}.{char_index}"
            current_index += len(line) + 1
        return "1.0"

    def detect_homopolymers(self):
        sequence = self.seqDict.get(self.seq_name)
        hp_list = SequenceAnalyzer.find_homopolymers(sequence)
        hps_print_format = self.hp_print(hp_list)
        hp_seq_formatted = self.convert_hp_to_lowercase(hp_list)
        self.display_full_sequence(hp_seq_formatted)
        self.display_hp(hps_print_format)

    def convert_hp_to_lowercase(self, hps=None):
        seq = self.seqDict.get(self.seq_name)
        if hps:
            for i in range(0, len(hps)):
                hp = hps[i]
                start = int(hp[1])
                end = int(hp[2])
                lowercase = seq[start:end].lower()
                if start == 0:
                    seq = lowercase + seq[end:]
                else:
                    seq = seq[0:start] + lowercase + seq[end:]
            formatted_seq = ""
            for i in range(0, len(seq), 80):
                formatted_seq += "  " + seq[i:i+80] + " \n"
        return formatted_seq

    def hp_print(self, hps):
        hp_print = "Homopolymers: \n"
        for i in range(0, len(hps)):
            hp = hps[i]
            type_hp = hp[0][0]
            start = int(hp[1])
            end = int(hp[2])
            length = end - start
            hp_print += str(i + 1) + "=" + type_hp + ":" + str(start) + "-" + str(end) + "_" + str(length) + "\n"
        return hp_print

    def clear_file(self):
        self.motif_text.delete(1.0, tk.END)
        self.cpg_text.delete(1.0, tk.END)
        self.hp_text.delete(1.0, tk.END)

if __name__ == "__main__":
    app = FastaViewer()
    app.mainloop()
