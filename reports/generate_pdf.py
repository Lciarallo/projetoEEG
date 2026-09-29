"""PDF local autocontido; erros de renderização encerram com código não zero."""
from pathlib import Path
import argparse
import base64
from datetime import date
import html
import mimetypes
import re
import shutil
import subprocess
import tempfile
import markdown

ROOT = Path(__file__).resolve().parents[1]
MD_PATH = ROOT/'reports'/'relatorio_cientifico.md'
PDF_OUT = ROOT/'reports'/'Relatorio_Cientifico_EEG_Binaural.pdf'


def render_html(md_content, md_path=MD_PATH):
    def image(match):
        alt, reference = match.groups()
        path = (md_path.parent/reference).resolve()
        if not path.is_file():
            raise FileNotFoundError(f'Figura ausente: {path}')
        mime = mimetypes.guess_type(path.name)[0] or 'application/octet-stream'
        encoded = base64.b64encode(path.read_bytes()).decode('ascii')
        return f'<figure><img src="data:{mime};base64,{encoded}" alt="{html.escape(alt)}"><figcaption>{html.escape(alt)}</figcaption></figure>'
    content = re.sub(r'!\[(.*?)\]\((.*?)\)',image,md_content)
    body = markdown.markdown(content,extensions=['tables','fenced_code'])
    reference_heading = '<h2>Referências verificadas</h2>'
    if reference_heading in body:
        body = body.replace(reference_heading, '<section class="references">'+reference_heading)+'</section>'
    return '''<!doctype html><html lang="pt-BR"><meta charset="utf-8">
<title>Áudio binaural e relaxamento: relatório científico</title><style>
@page {size:A4; margin:15mm 14mm 15mm; @bottom-right {content:"Página " counter(page);font:9pt sans-serif;color:#536172;}}
body {font:9.7pt/1.4 "DejaVu Sans",sans-serif;color:#202c38;}
h1 {font-size:21pt;line-height:1.2;color:#163b56;}
h2 {font-size:13pt;color:#163b56;margin-top:16pt;border-bottom:1px solid #ccd7df;padding-bottom:5pt;break-after:avoid;}
p {orphans:3;widows:3;}
table {width:100%;border-collapse:collapse;font-size:8pt;margin:9pt 0;break-inside:avoid;}
thead {display:table-header-group;} tr {break-inside:avoid;}
th,td {padding:4pt;border:1px solid #ccd7df;text-align:left;}
th {background:#e9f0f5;} figure {margin:12pt 0;break-inside:avoid;}
img {display:block;max-width:100%;max-height:205mm;margin:auto;}
figcaption {font-size:8pt;color:#536172;margin-top:5pt;text-align:center;}
code {font:8pt monospace;overflow-wrap:anywhere;} a {color:#215b83;overflow-wrap:anywhere;}
pre {white-space:pre-wrap;overflow-wrap:anywhere;}
li {margin-bottom:3pt;} .stamp {font-size:8pt;color:#536172;} .references {break-inside:avoid;}
</style><body>'''+f'<div class="stamp">Projeto EEG | Gerado em {date.today():%d/%m/%Y}</div>'+body+'</body></html>'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--markdown', type=Path, default=MD_PATH)
    parser.add_argument('--output', type=Path, default=PDF_OUT)
    args = parser.parse_args()
    browser = shutil.which('chromium-browser') or shutil.which('chromium') or shutil.which('google-chrome')
    if browser is None:
        raise RuntimeError('Instale Chromium ou Google Chrome para gerar o PDF')
    document = render_html(args.markdown.read_text(encoding='utf-8'), args.markdown.resolve())
    # Saída temporária impede que um PDF antigo seja confundido com geração bem-sucedida.
    with tempfile.TemporaryDirectory(prefix='eeg-report-') as directory:
        directory = Path(directory)
        source = directory/'report.html'
        target = directory/'report.pdf'
        source.write_text(document,encoding='utf-8')
        result = subprocess.run([browser,'--headless','--disable-gpu',
                    f'--user-data-dir={directory / "profile"}', '--no-pdf-header-footer',
                    f'--print-to-pdf={target}',source.as_uri()],capture_output=True,text=True,timeout=90)
        if result.returncode or not target.is_file() or target.stat().st_size < 1000:
            raise RuntimeError(f'Falha ao gerar PDF: {result.stderr[-2000:]}')
        args.output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(target,args.output)
        if args.output.resolve() == PDF_OUT.resolve():
            shutil.copyfile(target,ROOT/PDF_OUT.name)
    print(f'PDF atualizado: {args.output} ({args.output.stat().st_size} bytes).')


if __name__ == '__main__':
    main()
