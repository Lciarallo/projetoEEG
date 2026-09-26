import os
import re
import base64
import markdown
import subprocess

MD_PATH = "/home/lciarallo/.gemini/antigravity/brain/49f57a52-891c-4c0f-9afe-6958a0e4ea3f/relatorio_cientifico_definitivo_alta_confiabilidade.md"
HTML_OUT = "/tmp/relatorio_eeg.html"
PDF_OUT = "/home/lciarallo/projetoEEG/Relatorio_Cientifico_EEG_Binaural.pdf"
IMAGES_DIR = "/home/lciarallo/.gemini/antigravity/brain/49f57a52-891c-4c0f-9afe-6958a0e4ea3f"

def encode_image(img_path):
    with open(img_path, "rb") as f:
        data = f.read()
    ext = os.path.splitext(img_path)[1].lower().replace(".", "")
    return f"data:image/{ext};base64,{base64.b64encode(data).decode('utf-8')}"

def main():
    with open(MD_PATH, "r", encoding="utf-8") as f:
        md_content = f.read()

    # Substituir caminhos de imagens locais por base64 embutido para renderização perfeita
    def replace_img(match):
        alt = match.group(1)
        path = match.group(2)
        # Se for caminho local ou no brain, pega do resultados
        basename = os.path.basename(path)
        local_path = os.path.join(IMAGES_DIR, basename)
        if os.path.exists(local_path):
            b64 = encode_image(local_path)
            return f'<div class="figure-container"><img src="{b64}" alt="{alt}"/><div class="figure-caption">{alt}</div></div>'
        return match.group(0)

    md_content = re.sub(r'!\[(.*?)\]\((.*?)\)', replace_img, md_content)

    # Converter Markdown para HTML
    html_body = markdown.markdown(md_content, extensions=['tables', 'fenced_code', 'nl2br'])

    # Template HTML com CSS estilizado para PDF acadêmico/científico
    html_doc = f"""<!DOCTYPE html>
<html lang="pt-BR">
<head>
<meta charset="UTF-8">
<title>Investigação de EEG sob Batimentos Binaurais: Estudo Piloto Exploratório</title>
<style>
    @page {{
        size: A4;
        margin: 20mm 15mm 20mm 15mm;
        @bottom-right {{
            content: "Página " counter(page);
            font-size: 9pt;
            color: #718096;
        }}
    }}
    body {{
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
        color: #1a202c;
        line-height: 1.6;
        font-size: 11pt;
        background: #ffffff;
        margin: 0;
        padding: 0;
    }}
    h1 {{
        color: #1a365d;
        font-size: 20pt;
        border-bottom: 2px solid #2b6cb0;
        padding-bottom: 8px;
        margin-top: 0;
        margin-bottom: 15px;
        page-break-after: avoid;
    }}
    h2 {{
        color: #2b6cb0;
        font-size: 14pt;
        border-bottom: 1px solid #e2e8f0;
        padding-bottom: 5px;
        margin-top: 25px;
        margin-bottom: 12px;
        page-break-after: avoid;
    }}
    h3 {{
        color: #2d3748;
        font-size: 12pt;
        margin-top: 18px;
        margin-bottom: 8px;
        page-break-after: avoid;
    }}
    p {{
        margin-top: 0;
        margin-bottom: 10px;
        text-align: justify;
    }}
    table {{
        width: 100%;
        border-collapse: collapse;
        margin: 15px 0;
        font-size: 9.5pt;
        page-break-inside: avoid;
    }}
    th, td {{
        border: 1px solid #cbd5e0;
        padding: 7px 10px;
        text-align: left;
    }}
    th {{
        background-color: #ebf8ff;
        color: #2b6cb0;
        font-weight: 600;
    }}
    tr:nth-child(even) {{
        background-color: #f7fafc;
    }}
    ul, ol {{
        margin-top: 0;
        margin-bottom: 12px;
        padding-left: 22px;
    }}
    li {{
        margin-bottom: 6px;
    }}
    .figure-container {{
        text-align: center;
        margin: 18px 0;
        page-break-inside: avoid;
    }}
    img {{
        max-width: 95%;
        height: auto;
        border-radius: 4px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        border: 1px solid #e2e8f0;
    }}
    .figure-caption {{
        font-size: 9pt;
        color: #4a5568;
        margin-top: 6px;
        font-style: italic;
    }}
    code {{
        background-color: #edf2f7;
        color: #805ad5;
        padding: 2px 5px;
        border-radius: 3px;
        font-family: monospace;
        font-size: 9.5pt;
    }}
    .header-box {{
        background: linear-gradient(135deg, #ebf8ff 0%, #edf2f7 100%);
        padding: 12px 18px;
        border-radius: 6px;
        border-left: 4px solid #3182ce;
        margin-bottom: 20px;
    }}
    .header-box p {{
        margin: 0;
        font-size: 10pt;
        color: #2d3748;
    }}
    hr {{
        border: 0;
        border-top: 1px solid #e2e8f0;
        margin: 20px 0;
    }}
</style>
</head>
<body>
<div class="header-box">
    <p><strong>Projeto:</strong> Análise de EEG e Neuromodulação Acústica (OpenBCI Ganglion - 4 Canais)</p>
    <p><strong>Referência:</strong> Reanálise e Expansão Metodológica de Gao et al. (2014, <em>Int. J. Psychophysiol.</em>)</p>
    <p><strong>Data de Geração:</strong> 25/09/2026</p>
</div>
{html_body}
</body>
</html>
"""

    with open(HTML_OUT, "w", encoding="utf-8") as f:
        f.write(html_doc)
    print("HTML gerado em:", HTML_OUT)

    # Executar Chromium em modo headless para gerar PDF
    cmd = [
        "chromium-browser",
        "--headless",
        "--disable-gpu",
        "--no-sandbox",
        f"--print-to-pdf={PDF_OUT}",
        HTML_OUT
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode == 0:
        size_kb = os.path.getsize(PDF_OUT) / 1024
        print(f"PDF gerado com sucesso em: {PDF_OUT} ({size_kb:.1f} KB)")
    else:
        print("Erro ao gerar PDF:", res.stderr)

if __name__ == "__main__":
    main()
