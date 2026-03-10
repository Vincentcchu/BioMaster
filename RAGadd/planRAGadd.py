import sys
import os
import yaml

# Add the parent directory to Python path 
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from agents.Knowledge import generate_workflow_markdown, convert_markdown_to_json

if __name__ == "__main__":
    # Load API key from config.yaml
    config_path = "/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/agents/biomaster/config.yaml"
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    my_openai_key = config['api']['main']['key']
    my_openai_url = config['api']['main']['base_url']
    
    # my_content_path = "https://nf-co.re/oncoanalyser/1.0.0/"  # 可以是PDF路径或网页URL
    my_content_path = "/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/agents/biomaster/RAGadd/s13073-017-0467-4.pdf"  # 可以是PDF路径或网页URL
    # my_content_path = "/cs/student/projects2/aisd/2024/shekchu/projects/cell_agents/biomaster/RAGadd/bbaf207.pdf"  # 可以是PDF路径或网页URL
    
    # 生成Markdown格式的工作流
    markdown_file = "plan_knowledge.md"
    generate_workflow_markdown(my_content_path, my_openai_key, my_openai_url, markdown_file)
    
    # 将Markdown转换为JSON并保存
    json_file = "./doc/Plan_Knowledge.json"
    convert_markdown_to_json(markdown_file, json_file)
    
    print("工作流程提取和转换完成！")