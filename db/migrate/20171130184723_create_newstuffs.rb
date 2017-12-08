class CreateNewstuffs < ActiveRecord::Migration[5.1]
  def change
    create_table :newstuffs do |t|
	  t.text :description, null: false
	  t.timestamps null: false
    end
  end
end
