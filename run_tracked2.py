from agents.Biomaster import Biomaster
from langchain_core.messages import HumanMessage
from langchain_community.callbacks import get_openai_callback
import yaml
import json
import os
import sys
import time
from datetime import datetime

# Load configuration from YAML file
def load_config(config_path='config.yaml'):
    """Load YAML configuration file"""
    with open(config_path, 'r', encoding='utf-8') as file:
        return yaml.safe_load(file)

if __name__ == "__main__":
    # Parse command line arguments to get config file path
    config_path = sys.argv[1] if len(sys.argv) > 1 else 'config.yaml'
    
    # If config file doesn't exist, report error and exit
    if not os.path.exists(config_path):
        print(f"Error: Configuration file '{config_path}' does not exist")
        sys.exit(1)
    
    # Load configuration
    config = load_config(config_path)
    
    # Get API settings from configuration
    api_key = config['api']['main']['key']
    base_url = config['api']['main']['base_url']
    embedding_api_key = config['api']['embedding']['key']
    embedding_base_url = config['api']['embedding']['base_url']
    
    # Ollama settings
    use_ollama = config.get('biomaster', {}).get('use_ollama', False)
    ollama_config = config.get('api', {}).get('ollama', {})
    ollama_enabled = ollama_config.get('enabled', False)
    use_ollama = use_ollama and ollama_enabled
    ollama_base_url = ollama_config.get('base_url', 'http://localhost:11434')
    
    # Get model settings from configuration
    Model = config['models']['main']
    tool_model = config['models'].get('tool', Model)
    embedding_model = config['models']['embedding']
    
    # Get Biomaster settings from configuration
    excutor = config['biomaster']['executor']
    ids = config['biomaster']['id']
    generate_plan = config['biomaster'].get('generate_plan', True)
    
    # Get data and goal from configuration
    datalist = config['data']['files']
    goal = config['data']['goal']

    print(f"Using ollama: {use_ollama}")
    print(f"Ollama base URL: {ollama_base_url}")
    print(f"Model: {Model}")
    print(f"Embedding model: {embedding_model}")
    
    start_time = datetime.now()
    print(f"\n{'='*70}")
    print(f"RUN STARTED: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Run ID: {ids}")
    print(f"{'='*70}\n")

    try:
        # Wrap execution in token tracking callback
        with get_openai_callback() as cb:
            # Initialize Biomaster instance
            manager = Biomaster(
                api_key, 
                base_url,
                excutor=excutor,
                id=ids,
                Model=Model,
                embedding_model=embedding_model,
                tool_model=tool_model,
                embedding_base_url=embedding_base_url,
                embedding_api_key=embedding_api_key,
                use_ollama=use_ollama,
                ollama_base_url=ollama_base_url
            )
            
            # Execute PLAN and TASK
            if generate_plan:
                manager.execute_PLAN(goal, datalist)
                print("**********************************************************")

            PLAN_results_dict = manager.execute_TASK(datalist)
            print(PLAN_results_dict)
        
        # After execution, cb contains all token usage
        end_time = datetime.now()
        duration = end_time - start_time
        
        # Calculate cost manually based on model pricing
        # Update these rates based on your actual OpenAI pricing
        # https://openai.com/api/pricing/
        PRICING = {
            "gpt-4o": {"input": 2.50 / 1_000_000, "output": 10.00 / 1_000_000},
            "gpt-4o-mini": {"input": 0.150 / 1_000_000, "output": 0.600 / 1_000_000},
            "gpt-4-turbo": {"input": 10.00 / 1_000_000, "output": 30.00 / 1_000_000},
            "gpt-3.5-turbo": {"input": 0.50 / 1_000_000, "output": 1.50 / 1_000_000},
            "o1": {"input": 15.00 / 1_000_000, "output": 60.00 / 1_000_000},
            "o1-mini": {"input": 3.00 / 1_000_000, "output": 12.00 / 1_000_000},
            "o3-mini-2025-01-31": {"input": 1.10 / 1_000_000, "output": 4.40 / 1_000_000},
            "o3-mini": {"input": 1.10 / 1_000_000, "output": 4.40 / 1_000_000},
        }
        
        # Calculate cost based on model
        model_lower = Model.lower()
        manual_cost = 0.0
        cost_note = "(estimated, verify with OpenAI console)"
        
        for model_key, rates in PRICING.items():
            if model_key in model_lower:
                manual_cost = (cb.prompt_tokens * rates["input"]) + (cb.completion_tokens * rates["output"])
                cost_note = f"(based on {model_key} pricing)"
                break
        
        if manual_cost == 0.0:
            cost_note = "(unknown model, update PRICING dict)"
        
        # Create summary
        summary = {
            "run_id": ids,
            "config_file": config_path,
            "model": Model,
            "start_time": start_time.isoformat(),
            "end_time": end_time.isoformat(),
            "duration_seconds": duration.total_seconds(),
            "duration_formatted": str(duration),
            "token_usage": {
                "total_tokens": cb.total_tokens,
                "prompt_tokens": cb.prompt_tokens,
                "completion_tokens": cb.completion_tokens,
                "callback_reported_cost": cb.total_cost,  # What LangChain reports
                "calculated_cost": manual_cost,  # Manual calculation
                "cost_note": cost_note
            },
            "api_calls": cb.successful_requests
        }
        
        # Save summary
        os.makedirs('./output', exist_ok=True)
        summary_file = f'./output/{ids}_RUN_SUMMARY.json'
        with open(summary_file, 'w', encoding='utf-8') as f:
            json.dump(summary, f, indent=4)
        
        # Print summary
        print(f"\n{'='*70}")
        print(f"RUN COMPLETED: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"{'='*70}")
        print(f"Run ID:              {ids}")
        print(f"Model:               {Model}")
        print(f"Duration:            {duration}")
        print(f"\nToken Usage:")
        print(f"  Prompt tokens:     {cb.prompt_tokens:,}")
        print(f"  Completion tokens: {cb.completion_tokens:,}")
        print(f"  Total tokens:      {cb.total_tokens:,}")
        print(f"\nCost Estimates:")
        print(f"  Calculated cost:   ${manual_cost:.4f} {cost_note}")
        print(f"  Callback reported: ${cb.total_cost:.4f} (may be inaccurate)")
        print(f"  API calls:         {cb.successful_requests}")
        print(f"\n⚠️  Verify actual cost in OpenAI console")
        print(f"\nSummary saved to: {summary_file}")
        print(f"{'='*70}\n")
        
    except Exception as e:
        print(f"\nERROR during execution: {e}")
        raise